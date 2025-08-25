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

LMDB_PREFIX = os.path.join(settings.GWAS_LMDB_VAR_PVAL_DIR, "lmdb")  # shards: lmdb_chr{chrom}.mdb OR lmdb_chr{chrom}/

# ---------- metadata loaders (cached) ----------

_META_CACHE = {
    "traits": None,
    "n_traits": None,
    "metrics": None,
    "block_size": None,
    "trait_meta": None,  # phenocode -> dict(id, trait_group, trait_label)
}

def _open_env_any():
    """Open any shard to read global meta. Prefer chr1 if available."""
    # file-mode first
    candidate = f"{LMDB_PREFIX}_chr1.mdb"
    if os.path.exists(candidate):
        return lmdb.open(candidate, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=False)
    # dir-mode fallback
    candidate = f"{LMDB_PREFIX}_chr1"
    if os.path.isdir(candidate):
        return lmdb.open(candidate, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=True)
    # else pick any shard that exists
    parent = os.path.dirname(LMDB_PREFIX) or "."
    for name in os.listdir(parent):
        if name.startswith(os.path.basename(LMDB_PREFIX) + "_chr"):
            path = os.path.join(parent, name)
            if os.path.isfile(path) and path.endswith(".mdb"):
                return lmdb.open(path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=False)
            if os.path.isdir(path):
                return lmdb.open(path, readonly=True, lock=False, readahead=True, max_dbs=1, subdir=True)
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
            metrics_raw = txn.get(b"__meta__:metrics_json")
            block_raw   = txn.get(b"__meta__:block_size")
            metrics = json.loads(metrics_raw.decode()) if metrics_raw else ["p"]  # legacy
            block_size = int(block_raw.decode()) if block_raw else None
        _META_CACHE["traits"] = traits
        _META_CACHE["n_traits"] = n_traits
        _META_CACHE["metrics"] = metrics
        _META_CACHE["block_size"] = block_size
    finally:
        env.close()

    # Load phenotype metadata once
    pheno_csv = config("PHENO_FILE")
    df = pd.read_csv(pheno_csv)
    # be permissive with column names
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

# ---------- key helpers ----------

def _vkey_legacy(chrom, pos, ref, alt):
    return f"{chrom}:{pos}:{ref}:{alt}".encode()

def _vkey_block(chrom, pos, ref, alt, block_id: int):
    # block layout: "{chrom}:{pos}:{ref}:{alt}|{block_id}"
    return f"{chrom}:{pos}:{ref}:{alt}|{block_id}".encode()

# ---------- main function ----------

def extract_phenotype_results_for_variant(variant_id):
    """
    Returns:
    {
      "data": [ {id, trait_group, trait_label, log_pvalue, variant, chromosome, position, build, ref_allele}, ... ],
      "lastPage": null,
      "meta": { "build": ["GRCh37"] }
    }
    """
    _load_global_meta()
    traits = _META_CACHE["traits"]
    n_traits = _META_CACHE["n_traits"]
    metrics = _META_CACHE["metrics"]
    block_size = _META_CACHE["block_size"]
    trait_meta = _META_CACHE["trait_meta"]

    chrom, pos, ref, alt = convert_variant_id(variant_id)
    env, _ = _open_env_for_chrom(chrom)
    if env is None:
        return JsonResponse({"error": f"No LMDB shard for chrom {chrom}."}, status=404)

    try:
        # Fetch data either as blocks (new) or legacy single payload
        with env.begin() as txn:
            if block_size:  # block layout
                num_blocks = int(ceil(n_traits / float(block_size)))
                mats = []
                found_any = False
                for b in range(num_blocks):
                    key = _vkey_block(chrom, pos, ref, alt, b)
                    payload = txn.get(key)
                    if payload is None:
                        # empty block
                        mats.append(np.full((len(metrics), block_size), np.nan, dtype=np.float32))
                        continue
                    buf = ZSTD_DECOMP.decompress(payload) if USE_ZSTD else payload
                    arr = np.frombuffer(buf, dtype=np.float32).copy()
                    # Expected shape (M, BLOCK)
                    m = len(metrics)
                    if arr.size != m * block_size:
                        # tolerate partial last block
                        need = m * block_size
                        if arr.size < need:
                            pad = np.full(need - arr.size, np.nan, dtype=np.float32)
                            arr = np.concatenate([arr, pad], axis=0)
                        else:
                            arr = arr[:need]
                    mats.append(arr.reshape(len(metrics), block_size))
                    found_any = True
                if not found_any:
                    return JsonResponse({"error": f"Variant not found {chrom}:{pos}:{ref}:{alt}."}, status=404)
                mat_full = np.concatenate(mats, axis=1)[:, :n_traits]  # shape (M, n_traits)
            else:
                # legacy: single payload per variant; could be 1×N (p-only) or M×N
                key = _vkey_legacy(chrom, pos, ref, alt)
                payload = txn.get(key)
                if payload is None:
                    return JsonResponse({"error": f"Variant not found {chrom}:{pos}:{ref}:{alt}."}, status=404)
                buf = ZSTD_DECOMP.decompress(payload) if USE_ZSTD else payload
                arr = np.frombuffer(buf, dtype=np.float32).copy()
                m = len(metrics)
                if arr.size == n_traits:
                    # p-only
                    mat_full = np.full((m, n_traits), np.nan, dtype=np.float32)
                    mat_full[metrics.index("p"), :] = arr
                else:
                    mat_full = arr.reshape(m, n_traits)

        # Extract per-metric rows
        # Guarantee indices even if metric order changes
        def _get_row(name):
            return mat_full[metrics.index(name), :] if name in metrics else np.full(n_traits, np.nan, dtype=np.float32)

        p     = _get_row("p")
        beta  = _get_row("beta")
        se    = _get_row("se")
        af    = _get_row("af")

        # Build response list
        build = config("VITE_GENOME_BUILD", default="GRCh37")

        data_out = []
        for trait_code, pv, b, s, f in zip(traits, p, beta, se, af):
            if not np.isfinite(pv):
                continue
            meta = trait_meta.get(trait_code, {"id": None, "trait_group": None, "trait_label": trait_code})
            item = {
                "id": meta.get("id"),
                "trait_group": meta.get("trait_group"),
                "trait_label": meta.get("trait_label", trait_code),
                "log_pvalue": float(-np.log10(pv)) if pv > 0 else None,
                "beta": (float(b) if np.isfinite(b) else None),
                "se": (float(s) if np.isfinite(s) else None),
                "af": (float(f) if np.isfinite(f) else None)
            }
            data_out.append(item)

        data_out.sort(
            key=lambda x: (
                "" if x.get("trait_group") is None else str(x["trait_group"]),
                -(x.get("log_pvalue") if x.get("log_pvalue") is not None else float("-inf")),
            )
        )

        return {
            "data": data_out,
            "lastPage": None,
            "meta": {"build": [build]}
        }
    finally:
        env.close()
