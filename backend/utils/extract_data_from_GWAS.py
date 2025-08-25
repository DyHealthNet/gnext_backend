from backend.utils.converters import convert_variant_id
import os
import json
import lmdb
import numpy as np
import pandas as pd
from math import ceil
from django.conf import settings
from django.http import JsonResponse
from decouple import config
import logging
import contextlib
import struct, msgpack
import time
import re
import pysam

logger = logging.getLogger("backend")

# Optional zstd
try:
    import zstandard as zstd
    USE_ZSTD = True
    ZSTD_DECOMP = zstd.ZstdDecompressor()
except Exception:
    USE_ZSTD = False
    ZSTD_DECOMP = None

# Match the builder output path and naming
LMDB_PREFIX = os.path.join("test", "lmdb_p_f32")  # shards: lmdb_p_f32_chr{chrom}(.mdb or dir)

# ---------- metadata cache ----------

_META = {
    "traits": None,
    "n_traits": None,
    "trait_meta": None,  # phenocode -> {id, trait_group, trait_label}
}

def _open_env_any():
    """Open any shard to read global meta. Prefer chr1 if available."""
    # dir mode
    path = f"{LMDB_PREFIX}_chr1"
    return lmdb.open(path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=True), path

def _load_global_meta():
    env, path = _open_env_any()
    if env is None:
        raise RuntimeError("No LMDB shards found to read metadata.")
    logger.info(f"[GWAS LMDB] reading meta from: {path}")

    with env.begin() as txn:
        for key, val in txn.cursor():
            if key.startswith(b"__meta__"):
                try:
                    decoded_val = val.decode()
                except Exception:
                    decoded_val = val
                logger.info(f"{key!r}: {decoded_val}")

    env.close()

    try:
        with env.begin() as txn:
            traits = json.loads(txn.get(b"__meta__:traits_json").decode())
            n_traits = int(txn.get(b"__meta__:n_traits").decode())
        _META["traits"] = traits
        _META["n_traits"] = n_traits
    finally:
        env.close()

    # phenotype metadata (labels, groups) from PHENO_FILE
    pheno_csv = config("PHENO_FILE")
    df = pd.read_csv(pheno_csv)
    col_id = "id" if "id" in df.columns else "phenocode"
    col_group = "trait_group" if "trait_group" in df.columns else ("category" if "category" in df.columns else None)
    col_label = "trait_label" if "trait_label" in df.columns else ("description" if "description" in df.columns else None)

    trait_meta = {}
    for _, row in df.iterrows():
        code = str(row["phenocode"])
        meta = {"id": None, "trait_group": None, "trait_label": code}
        # id
        try:
            meta["id"] = int(row[col_id])
        except Exception:
            meta["id"] = row[col_id] if col_id in df.columns else code
        # group
        if col_group and pd.notna(row.get(col_group, None)):
            meta["trait_group"] = str(row[col_group])
        # label
        if col_label and pd.notna(row.get(col_label, None)):
            meta["trait_label"] = str(row[col_label])
        trait_meta[code] = meta
    _META["trait_meta"] = trait_meta

# ---------- env open helpers ----------

_ENV_CACHE = {}  # chrom_path -> lmdb.Environment

def _open_env_for_chrom(chrom: str):
    """Return a cached read-only env for the chromosome, or (None, None) if missing."""
    file_path = f"{LMDB_PREFIX}_chr{chrom}.mdb"
    dir_path = f"{LMDB_PREFIX}_chr{chrom}"

    key = None
    if os.path.exists(file_path):
        key = file_path
        subdir = False
    elif os.path.isdir(dir_path):
        key = dir_path
        subdir = True
    else:
        return None, None

    env = _ENV_CACHE.get(key)
    if env is None:
        env = lmdb.open(key, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=subdir)
        _ENV_CACHE[key] = env
    return env, subdir

# ---------- keys (p-only layout) ----------

def _vkey(chrom, pos, ref, alt):
    # single payload per variant: float32 vector length n_traits (p-only)
    return f"{chrom}:{pos}:{ref}:{alt}".encode()

# ---------- main function (p-only) ----------

def extract_phenotype_results_for_variant(variant_id):
    """
    Return p-only results across all traits for a given variant.
    Response:
    {
      "data": [ {id, trait_group, trait_label, log_pvalue, beta=None, se=None, af=None}, ... ],
      "lastPage": null,
      "meta": { "build": ["GRCh37"] }
    }
    """
    _load_global_meta()
    traits = _META["traits"]
    n_traits = _META["n_traits"]
    trait_meta = _META["trait_meta"]

    chrom, pos, ref, alt = convert_variant_id(variant_id)
    env, _ = _open_env_for_chrom(chrom)
    if env is None:
        return JsonResponse({"error": f"No LMDB shard for chrom {chrom}."}, status=404)

    # fetch p vector
    with env.begin() as txn:
        payload = txn.get(_vkey(chrom, pos, ref, alt))
        if payload is None:
            return JsonResponse({"error": f"Variant not found {chrom}:{pos}:{ref}:{alt}."}, status=404)
        buf = ZSTD_DECOMP.decompress(payload) if USE_ZSTD else payload
        arr = np.frombuffer(buf, dtype=np.float32)
        # tolerate minor mismatches
        if arr.size < n_traits:
            pad = np.full(n_traits - arr.size, np.nan, dtype=np.float32)
            p_vec = np.concatenate([arr, pad], axis=0)
        elif arr.size > n_traits:
            p_vec = arr[:n_traits]
        else:
            p_vec = arr

    build = config("VITE_GENOME_BUILD", default="GRCh37")

    # build output
    data_out = []
    for trait_code, pv in zip(traits, p_vec):
        if not np.isfinite(pv):
            continue
        meta = trait_meta.get(trait_code, {"id": None, "trait_group": None, "trait_label": trait_code})
        data_out.append({
            "id": meta.get("id"),
            "trait_group": meta.get("trait_group"),
            "trait_label": meta.get("trait_label", trait_code),
            "log_pvalue": pv,
            "beta": None,
            "se": None,
            "af": None,
        })

    # sort by group then by significance
    data_out.sort(
        key=lambda x: (
            "" if x.get("trait_group") is None else str(x["trait_group"]),
            -(x.get("log_pvalue") if x.get("log_pvalue") is not None else float("-inf")),
        )
    )

    return {
        "data": data_out,
        "lastPage": None,
        "meta": {"build": [build]},
    }


def extract_variants_for_range(filename, chr, start, end, pval_cutoff=1.0, max_rows=10000):
    # Build path to gzipped/tabix-indexed GWAS file
    norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(filename))
    norm_filepath = os.path.join(settings.GWAS_NORM_DIR, norm_filename + ".gz")

    # Initialize containers
    rows = []
    location = []
    ref = []
    alt = []
    external_ids = []
    allele_frequencies = []

    lmdb_path = os.path.join(settings.GWAS_NORM_DIR, f"lmdb_sorted_{config('VITE_GENOME_BUILD')}") + "/data.mdb"
    db_handles = {}

    try:
        lmdb_env = lmdb.open(
            lmdb_path,
            map_size=1024 ** 4,
            max_dbs=25,
            readonly=True,
            lock=False,
            readahead=True,
            subdir=False
        )
        logger.debug("LMDB environment opened")
    except Exception as e:
        logger.warning(f"LMDB not available, RSIDs will not be updated: {e}")
        lmdb_env = None

    tabix_file = pysam.TabixFile(norm_filepath)
    columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue',
               'pvalue', 'beta', 'stderr_beta', 'alt_allele_freq']

    with (lmdb_env.begin(buffers=True) if lmdb_env else contextlib.nullcontext()) as txn:
        try:
            for row in tabix_file.fetch(chr, start - 1, end):
                row = row.split("\t")
                if len(row) != len(columns):
                    logger.warning(f"Length of rows ({len(row) + 1}) does not match length of columns ({len(columns)}). Skipping malformed row: {row}")
                    continue

                row_dict = dict(zip(columns, row))

                # Convert numeric fields safely
                try:
                    pos = int(row_dict["pos"])
                    pvalue = float(row_dict["pvalue"]) if row_dict["pvalue"] != "." else None
                except ValueError:
                    logger.warning(f"Skipping row with invalid numeric values: {row}")
                    continue

                # Filter by p-value
                #if pvalue is not None and pvalue > pval_cutoff and pvalue != ".":
                if pvalue is not None and pvalue > pval_cutoff:
                    continue

                row_dict['variant_id'] = f"{chr}_{row_dict['pos']}_{row_dict['ref']}/{row_dict['alt']}"

                if lmdb_env and txn:
                    if chr not in db_handles:
                        db_handles[chr] = lmdb_env.open_db(chr.encode(), txn=txn)
                    db = db_handles[chr]

                    # Lookup RSID using LMDB
                    if db:
                        key_bytes = struct.pack('>I', pos)
                        value_bytes = txn.get(key_bytes, db=db)
                        if value_bytes:
                            refalt_to_rsid = msgpack.unpackb(value_bytes, raw=False)
                            refalt = f"{row_dict['ref']}/{row_dict['alt']}"
                            rsid_int = refalt_to_rsid.get(refalt)
                            if rsid_int is not None:
                                row_dict['rsid'] = f"rs{rsid_int}"

                rows.append(row_dict)
                location.append(int(row_dict["pos"]))
                ref.append(row_dict["ref"])
                alt.append(row_dict["alt"])
                external_ids.append(row_dict["rsid"])
                allele_frequencies.append(
                    float(row_dict["alt_allele_freq"]) if row_dict["alt_allele_freq"] != "." else None)

            temp_rows = len(rows)
            start_time = time.time()

            if max_rows is not None and len(rows) > max_rows:
                # Sort by pvalue ascending (None = worst)
                def sort_key(r):
                    try:
                        return float(r["pvalue"])
                    except (TypeError, ValueError):
                        return float("inf")

                sorted_indices = sorted(range(len(rows)), key=lambda i: sort_key(rows[i]))

                rows = [rows[i] for i in sorted_indices[:max_rows]]
                location = [location[i] for i in sorted_indices[:max_rows]]
                ref = [ref[i] for i in sorted_indices[:max_rows]]
                alt = [alt[i] for i in sorted_indices[:max_rows]]
                external_ids = [external_ids[i] for i in sorted_indices[:max_rows]]
                allele_frequencies = [allele_frequencies[i] for i in sorted_indices[:max_rows]]

            end_time = time.time()
            elapsed_time = end_time - start_time
            logger.debug(f"FINISHED FILTERING in {elapsed_time:.2f} seconds from {temp_rows} rows")

            data = {
                "header": ['variant_id'] + columns,
                "rows": rows,
                "location": location,
                "ref": ref,
                "alt": alt,
                "external_ids": external_ids,
                "allele_frequencies": allele_frequencies,
            }

            return data

        except Exception as e:
            logger.error(f"Error fetching data for {chr}:{start}-{end}: {e}")
            return {"error": str(e)}