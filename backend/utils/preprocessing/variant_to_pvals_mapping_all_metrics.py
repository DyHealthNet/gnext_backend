import os
import re
import json
import time
import logging
import multiprocessing as mp
from collections import defaultdict

import pandas as pd
import numpy as np
import lmdb
import pysam
from decouple import config
from django.conf import settings

logger = logging.getLogger("backend")

# =========================
# Config and constants
# =========================

METRICS = ["p", "beta", "se", "af"]
M = len(METRICS)

COL = {
    "p": 5,
    "beta": 7,
    "se": 8,
    "af": 9,
}

MAX_WORKERS = int(config("MAX_WORKERS"))

# Batch: how many variants to buffer before flushing to LMDB
VARIANT_BATCH_SIZE = 20_000   # tune this depending on your RAM


# =========================
# Helpers
# =========================

def vkey(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}:{alt}".encode()


def parse_rec(line):
    f = line.rstrip().split("\t")
    try:
        pos = int(f[1])
    except Exception:
        return None
    chrom, ref, alt = f[0], f[3], f[4]

    def as_float(x):
        try:
            return float(x) if x != "." else np.nan
        except Exception:
            return np.nan

    vals = [
        as_float(f[COL["p"]])   if len(f) > COL["p"]   else np.nan,
        as_float(f[COL["beta"]])if len(f) > COL["beta"]else np.nan,
        as_float(f[COL["se"]])  if len(f) > COL["se"]  else np.nan,
        as_float(f[COL["af"]])  if len(f) > COL["af"]  else np.nan,
    ]
    return chrom, pos, ref, alt, vals


def init_env(lmdb_path, traits):
    env = lmdb.open(
        lmdb_path,
        map_size=1 << 40,
        subdir=True,
        writemap=False,
        metasync=False,
        sync=False,
        max_dbs=1,
        readahead=False,   # speeds up bulk writes
    )
    with env.begin(write=True) as txn:
        txn.put(b"__meta__:n_traits", str(len(traits)).encode())
        txn.put(b"__meta__:traits_json", json.dumps(traits).encode())
        txn.put(b"__meta__:metrics_json", json.dumps(METRICS).encode())
    return env


def put_mat(txn, key, mat):
    txn.put(key, mat.tobytes())


# =========================
# Worker: build one chromosome shard
# =========================

def build_chrom(args):
    chrom, lmdb_root, traits, t2p = args
    n_traits = len(traits)
    trait_index = {t: i for i, t in enumerate(traits)}
    lmdb_path = f"{lmdb_root}_chr{chrom}"
    env = init_env(lmdb_path, traits)

    txn = env.begin(write=True)
    variant_store = defaultdict(lambda: np.full((M, n_traits), np.nan, dtype=np.float32))
    count = 0
    flushed = 0

    try:
        # Iterate all traits once, fill matrices in memory
        for trait in traits:
            j = trait_index[trait]
            tbx = pysam.TabixFile(t2p[trait])
            try:
                it = tbx.fetch(chrom)
            except ValueError:
                tbx.close()
                continue

            for line in it:
                rec = parse_rec(line)
                if rec is None:
                    continue
                c, pos, ref, alt, vals = rec
                key = vkey(c, pos, ref, alt)
                mat = variant_store[key]
                for m, v in enumerate(vals):
                    if not np.isnan(v) and np.isnan(mat[m, j]):
                        mat[m, j] = v

                count += 1
                # Flush if memory batch exceeded
                if count % VARIANT_BATCH_SIZE == 0:
                    for k, mtx in variant_store.items():
                        put_mat(txn, k, mtx)
                    txn.commit()
                    txn = env.begin(write=True)
                    flushed += len(variant_store)
                    variant_store.clear()

            tbx.close()

        # Final flush
        for k, mtx in variant_store.items():
            put_mat(txn, k, mtx)
        txn.commit()
        flushed += len(variant_store)
        logger.info(f"chr{chrom}: trait {trait} stored in LMDB")



    except Exception:
        txn.abort()
        raise
    finally:
        env.sync()
        env.close()

    return chrom


# =========================
# Orchestrator
# =========================

def build_snp_pval_map_lmdb(lmdb_root_prefix):
    pheno_file = config("PHENO_FILE")
    pheno_dt = pd.read_csv(pheno_file)
    traits = pheno_dt["phenocode"].astype(str).tolist()
    paths = [
        os.path.join(
            settings.GWAS_NORM_DIR,
            re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(i)) + ".gz"
        )
        for i in pheno_dt["filename"].astype(str).tolist()
    ]
    t2p = dict(zip(traits, paths))

    # Chromosomes from the first file
    first_file = paths[0]
    with pysam.TabixFile(first_file) as tbx:
        chroms = list(tbx.contigs)

    logger.info(f"Building LMDB shards for chromosomes: {chroms}")
    with mp.Pool(processes=MAX_WORKERS) as pool:
        for done in pool.imap_unordered(
            build_chrom, [(c, lmdb_root_prefix, traits, t2p) for c in chroms]
        ):
            logger.info(f"...Finished LMDB creation of chromosome {done}...")


def setup_variant_pvals_mapping_lmdb():
    lmdb_path = os.path.join(settings.GWAS_LMDB_VAR_PVAL_DIR, "lmdb")
    os.makedirs(settings.GWAS_LMDB_VAR_PVAL_DIR, exist_ok=True)
    start_time = time.time()
    build_snp_pval_map_lmdb(lmdb_path)
    elapsed = time.time() - start_time
    logger.debug(f"Time taken to produce LMDB mapping lib: {elapsed:.2f} seconds")
    return lmdb_path
