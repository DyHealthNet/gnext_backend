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

logger = logging.getLogger("backend")

# Optional zstd
try:
    import zstandard as zstd
    USE_ZSTD = True
    ZSTD_DECOMP = zstd.ZstdDecompressor()
except Exception:
    USE_ZSTD = False
    ZSTD_DECOMP = None

LMDB_PREFIX = os.path.join(settings.GWAS_LMDB_VAR_PVAL_DIR, "lmdb")  # shards: lmdb_chr{chrom}(.mdb|/)

# ---------- metadata cache ----------

_META_CACHE = {
    "traits": None,
    "n_traits": None,
    "block_size": None,
    "trait_meta": None,   # phenocode -> {id, trait_group, trait_label}
}

def _open_env_any():
    """Open any shard to read global meta. Prefer chr1 if available."""
    # file mode
    path = f"{LMDB_PREFIX}_chr1.mdb"
    if os.path.exists(path):
        return lmdb.open(path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=False)
    # dir mode
    path = f"{LMDB_PREFIX}_chr1"
    if os.path.isdir(path):
        return lmdb.open(path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=True)
    # fallback: first shard found
    parent = os.path.dirname(LMDB_PREFIX) or "."
    base = os.path.basename(LMDB_PREFIX) + "_chr"
    for name in os.listdir(parent):
        if not name.startswith(base):
            continue
        p = os.path.join(parent, name)
        if os.path.isfile(p) and p.endswith(".mdb"):
            return lmdb.open(p, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=False)
        if os.path.isdir(p):
            return lmdb.open(p, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=True)
    return None

def _load_global_meta():
    if _META_CACHE["traits"] is not None:
        return
    env = _open_env_any()
    if env is None:
        raise RuntimeError("No LMDB shards found to read metadata.")
    try:
        with env.begin() as txn:
            traits = json.loads(txn.get(b"__meta__:traits_json").decode())
            n_traits = int(txn.get(b"__meta__:n_traits").decode())
            block_raw = txn.get(b"__meta__:block_size")  # may be absent in legacy
            block_size = int(block_raw.decode()) if block_raw else None
        _META_CACHE["traits"] = traits
        _META_CACHE["n_traits"] = n_traits
        _META_CACHE["block_size"] = block_size
    finally:
        env.close()

    # phenotype metadata (labels, groups)
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
    _META_CACHE["trait_meta"] = trait_meta

def _open_env_for_chrom(chrom: str):
    """Return (env, subdir_flag) for given chromosome, or (None, None) if missing."""
    file_path = f"{LMDB_PREFIX}_chr{chrom}.mdb"
    if os.path.exists(file_path):
        return lmdb.open(file_path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=False), False
    dir_path = f"{LMDB_PREFIX}_chr{chrom}"
    if os.path.isdir(dir_path):
        return lmdb.open(dir_path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=True), True
    return None, None

# ---------- key helpers (p-only layout) ----------

def _vkey_legacy(chrom, pos, ref, alt):
    # single payload per variant: float32 vector length n_traits (p-only)
    return f"{chrom}:{pos}:{ref}:{alt}".encode()

def _vkey_block(chrom, pos, ref, alt, block_id: int):
    # block payload per variant: float32 vector length block_size (p-only)
    return f"{chrom}:{pos}:{ref}:{alt}|{block_id}".encode()

# ---------- main function (p-only) ----------

def extract_phenotype_results_for_variant(variant_id):
    """
    Returns p-only results across all traits for a given variant.
    Response:
    {
      "data": [ {id, trait_group, trait_label, log_pvalue, beta=None, se=None, af=None}, ... ],
      "lastPage": null,
      "meta": { "build": ["GRCh37"] }
    }
    """
    _load_global_meta()
    traits = _META_CACHE["traits"]
    n_traits = _META_CACHE["n_traits"]
    block_size = _META_CACHE["block_size"]
    trait_meta = _META_CACHE["trait_meta"]

    chrom, pos, ref, alt = convert_variant_id(variant_id)
    env, _ = _open_env_for_chrom(chrom)
    if env is None:
        return JsonResponse({"error": f"No LMDB shard for chrom {chrom}."}, status=404)

    try:
        with env.begin() as txn:
            if block_size:
                # Blocked layout: concatenate blocks into one 1D float32 vector of length n_traits
                num_blocks = int(ceil(n_traits / float(block_size)))
                segments = []
                found_any = False
                for b in range(num_blocks):
                    key = _vkey_block(chrom, pos, ref, alt, b)
                    payload = txn.get(key)
                    if payload is None:
                        # missing block segment → fill with NaNs
                        segments.append(np.full(block_size, np.nan, dtype=np.float32))
                        continue
                    buf = ZSTD_DECOMP.decompress(payload) if USE_ZSTD else payload
                    arr = np.frombuffer(buf, dtype=np.float32).copy()  # writable
                    if arr.size < block_size:
                        pad = np.full(block_size - arr.size, np.nan, dtype=np.float32)
                        arr = np.concatenate([arr, pad], axis=0)
                    elif arr.size > block_size:
                        arr = arr[:block_size]
                    segments.append(arr)
                    found_any = True
                if not found_any:
                    return JsonResponse({"error": f"Variant not found {chrom}:{pos}:{ref}:{alt}."}, status=404)
                p_vec = np.concatenate(segments, axis=0)[:n_traits]  # shape (n_traits,)
            else:
                # Legacy layout: single payload per variant is 1D float32 vector of length n_traits
                key = _vkey_legacy(chrom, pos, ref, alt)
                payload = txn.get(key)
                if payload is None:
                    return JsonResponse({"error": f"Variant not found {chrom}:{pos}:{ref}:{alt}."}, status=404)
                buf = ZSTD_DECOMP.decompress(payload) if USE_ZSTD else payload
                arr = np.frombuffer(buf, dtype=np.float32).copy()
                if arr.size != n_traits:
                    # tolerate minor mismatches by pad/trim
                    if arr.size < n_traits:
                        pad = np.full(n_traits - arr.size, np.nan, dtype=np.float32)
                        arr = np.concatenate([arr, pad], axis=0)
                    else:
                        arr = arr[:n_traits]
                p_vec = arr

        build = config("VITE_GENOME_BUILD", default="GRCh37")

        # Build output
        data_out = []
        for trait_code, pv in zip(traits, p_vec):
            if not np.isfinite(pv):
                continue
            meta = trait_meta.get(trait_code, {"id": None, "trait_group": None, "trait_label": trait_code})
            data_out.append({
                "id": meta.get("id"),
                "trait_group": meta.get("trait_group"),
                "trait_label": meta.get("trait_label", trait_code),
                "log_pvalue": float(-np.log10(pv)) if pv > 0 else None,
                "beta": None,
                "se": None,
                "af": None,
            })

        # Sort by group then by significance
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
    finally:
        env.close()