#!/usr/bin/env python3
import os, re, json, time, sys, logging, multiprocessing as mp
from collections import OrderedDict

import numpy as np
import pandas as pd
import lmdb
import pysam
from decouple import config
from django.conf import settings
import zstandard as zstd

# =========================
# Config
# =========================
PHENO_FILE    = config("PHENO_FILE")                       # CSV with columns: phenocode, filename
GWAS_NORM_DIR = settings.GWAS_NORM_DIR                     # Folder with .gz + .tbi per trait
OUT_DIR       = settings.GWAS_LMDB_VAR_PVAL_DIR            # Output directory
MAX_WORKERS   = int(config("MAX_WORKERS", default="8"))
P_COL         = 5                                          # 0-based p-value column

# LMDB settings
# 8M variants × 7k traits × 4 bytes ≈ 224 GB raw across all chromosomes before compression/overhead.
MAP_SIZE_BYTES = 512 * (1 << 30)                           # 512 GB
LMDB_FLAGS = dict(subdir=True, writemap=False, metasync=False,
                  sync=False, max_dbs=1, readahead=False, map_size=MAP_SIZE_BYTES)

# Compression
ZSTD_LEVEL  = 7
ZSTD_COMP   = zstd.ZstdCompressor(level=ZSTD_LEVEL)
ZSTD_DECOMP = zstd.ZstdDecompressor()

# Cache: each vector is n_traits*4 bytes raw; with Python/NumPy overhead ~60–80 KB at 7k traits
CACHE_CAPACITY    = int(config("CACHE_CAPACITY", default="20000"))   # good default for 64 GB with 8 workers
COMMIT_EVERY_PUTS = 200_000

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(processName)s %(levelname)s %(message)s",
    stream=sys.stdout,
    force=True,
)
logger = logging.getLogger("lmdb_p_f32")

# =========================
# Helpers
# =========================
def variant_key(chrom: str, pos: int, ref: str, alt: str) -> bytes:
    return f"{chrom}:{pos}:{ref}:{alt}".encode()

def parse_line(line: str):
    f = line.rstrip().split("\t")
    if len(f) < 5:
        return None
    try:
        pos = int(f[1])
    except Exception:
        return None
    return f[0], pos, f[3], f[4]

def get_pvalue(line: str) -> float:
    f = line.rstrip().split("\t")
    if P_COL >= len(f):
        return np.nan
    x = f[P_COL]
    if x == "." or x == "":
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan

def empty_vector(n_traits: int) -> np.ndarray:
    # float32 with NaN sentinel for missing
    return np.full((n_traits,), np.nan, dtype=np.float32)

def compress_vec(arr: np.ndarray) -> bytes:
    # Ensure C-contiguous
    b = np.ascontiguousarray(arr).view(np.uint8)
    return ZSTD_COMP.compress(b)

def decompress_vec(buf: bytes, n_traits: int) -> np.ndarray:
    # IMPORTANT: .copy() to make it writeable
    raw = ZSTD_DECOMP.decompress(buf, max_output_size=n_traits * 4)
    return np.frombuffer(raw, dtype=np.float32, count=n_traits).copy()

# =========================
# LRU write-back cache
# =========================
class LRUCache:
    def __init__(self, capacity: int, n_traits: int, env: lmdb.Environment):
        self.capacity = capacity
        self.n_traits = n_traits
        self.env = env
        self.store = OrderedDict()   # key -> (np.ndarray, dirty)
        self.txn = env.begin(write=True)
        self.put_count = 0

    def _commit_maybe(self):
        self.put_count += 1
        if self.put_count % COMMIT_EVERY_PUTS == 0:
            self.txn.commit()
            self.txn = self.env.begin(write=True)

    def _evict_one(self):
        k, (vec, dirty) = self.store.popitem(last=False)
        if dirty:
            self.txn.put(k, compress_vec(vec))
            self._commit_maybe()

    def get_for_update(self, key: bytes) -> np.ndarray:
        if key in self.store:
            vec, dirty = self.store.pop(key)
            self.store[key] = (vec, dirty)
            return vec
        if len(self.store) >= self.capacity:
            self._evict_one()
        buf = self.txn.get(key)
        vec = empty_vector(self.n_traits) if buf is None else decompress_vec(buf, self.n_traits)
        if not vec.flags.writeable:  # belt & suspenders
            vec = vec.copy()
        self.store[key] = (vec, False)
        return vec

    def mark_dirty(self, key: bytes):
        vec, _ = self.store.pop(key)
        self.store[key] = (vec, True)

    def flush_all(self):
        for k, (vec, dirty) in list(self.store.items()):
            if dirty:
                self.txn.put(k, compress_vec(vec))
        self.txn.commit()
        self.txn = None
        self.store.clear()

# =========================
# Build per chromosome
# =========================
def init_env(lmdb_path: str, traits: list[str]):
    os.makedirs(os.path.dirname(lmdb_path), exist_ok=True)
    env = lmdb.open(lmdb_path, **LMDB_FLAGS)
    with env.begin(write=True) as txn:
        txn.put(b"__meta__:n_traits", str(len(traits)).encode())
        txn.put(b"__meta__:traits_json", json.dumps(traits).encode())
        meta = {
            "encoding": {"metric": "p", "dtype": "float32", "nan_is_missing": True},
            "codec": {"type": "zstd", "level": ZSTD_LEVEL},
            "p_col": P_COL,
        }
        txn.put(b"__meta__:info", json.dumps(meta).encode())
    return env

def build_chrom(args):
    chrom, lmdb_root, traits, t2p, cache_capacity = args
    n_traits = len(traits)
    env = init_env(f"{lmdb_root}_chr{chrom}", traits)
    cache = LRUCache(cache_capacity, n_traits, env)

    chrom_t0 = time.time()
    total_lines = 0
    total_updates = 0
    env_closed = False

    try:
        for j, trait in enumerate(traits):
            trait_t0 = time.time()
            lines_processed = 0
            updates = 0
            last_log_t = trait_t0
            next_lines_log = 1_000_000

            logger.info(f"chr{chrom} trait={trait} start ({j+1}/{len(traits)})")

            tbx = pysam.TabixFile(t2p[trait])
            try:
                it = tbx.fetch(chrom)
            except ValueError:
                tbx.close()
                logger.info(f"chr{chrom} trait={trait}: no records on this contig")
                continue

            for line in it:
                head = parse_line(line)
                if not head:
                    continue
                c, pos, ref, alt = head
                key = variant_key(c, pos, ref, alt)

                p = get_pvalue(line)
                if not np.isnan(p):
                    vec = cache.get_for_update(key)
                    if np.isnan(vec[j]):
                        vec[j] = np.float32(p)   # store raw p-value
                        cache.mark_dirty(key)
                        updates += 1
                lines_processed += 1

                # heartbeat
                now = time.time()
                if (now - last_log_t) >= 5.0 or lines_processed >= next_lines_log:
                    elapsed = now - trait_t0
                    rate = lines_processed / max(elapsed, 1e-6)
                    logger.info(
                        f"chr{chrom} trait={trait} progress: "
                        f"lines={lines_processed:,} updates={updates:,} "
                        f"elapsed={elapsed:.1f}s rate={rate:,.0f}/s"
                    )
                    last_log_t = now
                    next_lines_log = ((lines_processed // 1_000_000) + 1) * 1_000_000

            tbx.close()

            dt_trait = time.time() - trait_t0
            total_lines += lines_processed
            total_updates += updates
            logger.info(
                f"chr{chrom} trait={trait} done: lines={lines_processed:,} updates={updates:,} "
                f"elapsed={dt_trait:.1f}s ({j+1}/{len(traits)})"
            )

        cache.flush_all()
    except Exception:
        try:
            if cache.txn is not None:
                cache.txn.abort()
        except Exception:
            pass
        try:
            env.close()
            env_closed = True
        except Exception:
            pass
        raise
    finally:
        if not env_closed:
            try:
                env.sync()
            except Exception:
                pass
            try:
                env.close()
            except Exception:
                pass

    dt = time.time() - chrom_t0
    logger.info(f"chr{chrom} summary: lines={total_lines:,} updates={total_updates:,} time={dt:.1f}s")
    return chrom

# =========================
# Orchestrator
# =========================
def build_lmdb_neglogp_f32():
    """
    Kept name for backward compatibility. Writes p-only (float32) shards:
    OUT_DIR/'lmdb_p_f32_chr{chrom}'.
    """
    pheno = pd.read_csv(PHENO_FILE)
    traits = pheno["phenocode"].astype(str).tolist()
    paths = [
        os.path.join(
            GWAS_NORM_DIR,
            re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(i)) + ".gz"
        )
        for i in pheno["filename"].astype(str).tolist()
    ]
    t2p = dict(zip(traits, paths))

    first_file = paths[0]
    with pysam.TabixFile(first_file) as tbx:
        chroms = list(tbx.contigs)

    os.makedirs(OUT_DIR, exist_ok=True)
    lmdb_root = os.path.join(OUT_DIR, "lmdb_p_f32")  # renamed to reflect p-only
    workers = min(MAX_WORKERS, len(chroms))

    logger.info(f"Build LMDB p-only float32 + zstd; P_COL={P_COL}; traits={len(traits)}; workers={workers}")
    args = [(c, lmdb_root, traits, t2p, CACHE_CAPACITY) for c in chroms]

    if workers > 1:
        with mp.Pool(processes=workers) as pool:
            for done in pool.imap_unordered(build_chrom, args):
                logger.info(f"...finished chr{done}")
    else:
        for a in args:
            build_chrom(a)

    logger.info("Done.")
    return lmdb_root