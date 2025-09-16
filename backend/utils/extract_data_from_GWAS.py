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
import math
import re
import pysam
import heapq
import itertools
import subprocess
from backend.utils.pheno_cache import get_pheno_df, get_pheno_map
from backend.utils.preprocessing.zorp.zorp import sniffers
from backend.utils.preprocessing.snp_to_rsid_mapping import add_rsID_with_lmdb

logger = logging.getLogger("backend")


def _to_float(s):
    if s in (".", "NA", "", None):
        return None
    try:
        return float(s)
    except Exception:
        return None

def extract_variant_metrics(chr, pos, ref, alt):
    pheno_file = config("PHENO_FILE")
    pheno_map = get_pheno_map(pheno_file)  # cached in memory

    target_vid = f"{chr}:{pos}:{ref.upper()}:{alt.upper()}"

    # Paths to metric files
    metric_files = {
        "neg_log_pvalue": f"{settings.GWAS_CHR_BGZ_DIR}/chr{chr}/neg_log_pvalue.tsv.bgz",
        "beta": f"{settings.GWAS_CHR_BGZ_DIR}/chr{chr}/beta.tsv.bgz",
        "stderr_beta": f"{settings.GWAS_CHR_BGZ_DIR}/chr{chr}/stderr_beta.tsv.bgz",
        "alt_allele_freq": f"{settings.GWAS_CHR_BGZ_DIR}/chr{chr}/alt_allele_freq.tsv.bgz",
    }

    # Get trait names from header
    with pysam.TabixFile(metric_files["neg_log_pvalue"]) as tbx_p:
        header_line = tbx_p.header[0].lstrip("#")
        trait_names = header_line.split("\t")[5:]
        cand_lines = list(tbx_p.fetch(chr, pos - 1, pos))
        if not cand_lines:
            return None
        # Find the line with matching vid
        p_line = next((ln for ln in cand_lines if ln.split("\t")[4] == target_vid), None)
        if p_line is None:
            return None

    # Fetch the matching vid from each metric file
    metric_data = {}
    for metric_name, path in metric_files.items():
        with pysam.TabixFile(path) as tbx:
            lines = list(tbx.fetch(chr, pos - 1, pos))
        match_line = next((ln for ln in lines if ln.split("\t")[4] == target_vid), None)
        if match_line:
            metric_data[metric_name] = match_line.split("\t")[5:]
        else:
            metric_data[metric_name] = ["." for _ in trait_names]

    results = []

    # compute min and max allele frequency
    af_arr = np.array([_to_float(v) for v in metric_data["alt_allele_freq"]], dtype=float)
    mask = ~np.isnan(af_arr)
    min_af = float(np.min(af_arr[mask])) if mask.any() else float("inf")
    max_af = float(np.max(af_arr[mask])) if mask.any() else float("-inf")

    for idx, trait_code in enumerate(trait_names):
        cat, desc = pheno_map.get(trait_code, ("", trait_code))
        results.append({
            "id": trait_code,
            "x": idx,
            "trait_group": cat,
            "trait_label": desc,
            "log_pvalue": _to_float(metric_data["neg_log_pvalue"][idx]),
            "pvalue": np.power(10, -_to_float(metric_data["neg_log_pvalue"][idx])),
            "beta": _to_float(metric_data["beta"][idx]),
            "stderr_beta": _to_float(metric_data["stderr_beta"][idx]),
            "alt_allele_freq":_to_float(metric_data["alt_allele_freq"][idx])
        })

    # results = [
    #     {**r, "x": i}
    #     for i, r in enumerate(sorted(
    #         results,
    #         key=lambda r: (
    #             r["trait_group"],
    #             -r["log_pvalue"] if r["log_pvalue"] is not None else float("-inf")
    #         )
    #     ))
    # ]
    return results, min_af, max_af

def extract_variants_for_range(filename, chr, start, end, pval_cutoff=1.0, max_rows=10000):
    # Build path to gzipped/tabix-indexed GWAS file
    norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(filename))
    norm_filepath = os.path.join(settings.GWAS_NORM_DIR, norm_filename + ".gz")

    lmdb_path = os.path.join(settings.GWAS_NORM_DIR, f"lmdb_sorted_{config('VITE_GENOME_BUILD')}") + "/data.mdb"
    db_handles = {}

    heap = []
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
    col_idx = {name: i for i, name in enumerate(columns)}
    idx_nlog = col_idx['neg_log_pvalue']
    neg_log_cutoff = -math.log10(pval_cutoff)

    with (lmdb_env.begin(buffers=True) if lmdb_env else contextlib.nullcontext()) as txn:
        try:
            for row in tabix_file.fetch(chr, start - 1, end):
                row = row.split("\t")
                if len(row) != len(columns):
                    logger.warning(
                        f"Length of rows ({len(row)}) does not match length of columns ({len(columns)}). Skipping malformed row: {row}")
                    continue

                neg_log_pval = float(row[idx_nlog]) if row[idx_nlog] != "." else None
                if neg_log_pval is None or neg_log_pval < neg_log_cutoff:
                    continue

                if len(heap) < max_rows:
                    heapq.heappush(heap, (neg_log_pval, row))  # directly store the neg_log_pval
                else:
                    if neg_log_pval > heap[0][0]:  # if current value is greater than the smallest in the heap
                        heapq.heapreplace(heap, (neg_log_pval, row))

            rows = [dict(zip(columns, r)) for _, r in heap]
            for r in rows:

                pos = r['pos']
                chr = r['chrom']
                ref = r['ref']
                alt = r['alt']

                r['variant_id'] = f"{chr}_{pos}_{ref}/{alt}"

                if lmdb_env and txn:
                    if chr not in db_handles:
                        db_handles[chr] = lmdb_env.open_db(chr.encode(), txn=txn)
                    db = db_handles[chr]

                    # Lookup RSID using LMDB
                    if db:
                        key_bytes = struct.pack('>I', int(pos))
                        value_bytes = txn.get(key_bytes, db=db)
                        if value_bytes:
                            refalt_to_rsid = msgpack.unpackb(value_bytes, raw=False)
                            refalt = f"{ref}/{alt}"
                            rsid_int = refalt_to_rsid.get(refalt)
                            if rsid_int is not None:
                                r['rsid'] = f"rs{rsid_int}"

            data = {
                "header": ['variant_id'] + columns,
                "rows": rows,
            }

            return data

        except Exception as e:
            logger.error(f"Error fetching data for {chr}:{start}-{end}: {e}")
            return {"error": str(e)}


def get_all_sign_variants_cutoff(filename, pval_cutoff=5e-8, max_rows=10000):
    # Normalize filename
    norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(filename))
    norm_filepath = os.path.join(settings.GWAS_NORM_DIR, norm_filename + ".gz")

    lmdb_path = os.path.join(settings.GWAS_NORM_DIR, f"lmdb_{config('VITE_GENOME_BUILD')}") + "/data.mdb"
    db_handles = {}

    heap = []

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
    except Exception:
        lmdb_env = None

    columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue',
               'pvalue', 'beta', 'stderr_beta', 'alt_allele_freq']
    col_idx = {name: i for i, name in enumerate(columns)}
    idx_nlog = col_idx['neg_log_pvalue']
    neg_log_cutoff = -math.log10(pval_cutoff)

    with (lmdb_env.begin(buffers=True) if lmdb_env else contextlib.nullcontext()) as txn:
        try:
            for row in stream_filtered_variants(norm_filepath, neg_log_cutoff=neg_log_cutoff):
                if len(row) != len(columns):
                    logger.warning(
                        f"Length of rows ({len(row)}) does not match length of columns ({len(columns)}). Skipping malformed row: {row}")
                    continue
                if row[0].startswith("#") or row[idx_nlog] == "neg_log_pvalue":
                    continue

                neg_log_pval = float(row[idx_nlog]) if row[idx_nlog] != "." else None
                if neg_log_pval is None or neg_log_pval < neg_log_cutoff:
                    continue

                if len(heap) < max_rows:
                    heapq.heappush(heap, (neg_log_pval, row))  # directly store the neg_log_pval
                else:
                    if neg_log_pval > heap[0][0]:  # if current value is greater than the smallest in the heap
                        heapq.heapreplace(heap, (neg_log_pval, row))

            rows = [dict(zip(columns, r)) for _, r in heap]
            for r in rows:

                pos = r['pos']
                chr = r['chrom']
                ref = r['ref']
                alt = r['alt']

                r['variant_id'] = f"{chr}_{pos}_{ref}/{alt}"

                if lmdb_env and txn:
                    if chr not in db_handles:
                        db_handles[chr] = lmdb_env.open_db(chr.encode(), txn=txn)
                    db = db_handles[chr]

                    # Lookup RSID using LMDB
                    if db:
                        key_bytes = struct.pack('I', int(pos))
                        #key_bytes = struct.pack('>I', pos)
                        value_bytes = txn.get(key_bytes, db=db)
                        if value_bytes:
                            refalt_to_rsid = msgpack.unpackb(value_bytes, raw=False)
                            refalt = f"{ref}/{alt}"
                            rsid_int = refalt_to_rsid.get(refalt)
                            if rsid_int is not None:
                                r['rsid'] = f"rs{rsid_int}"

            data = {
                "header": ['variant_id'] + columns,
                "rows": rows,
            }

            return data

        except Exception as e:
            logger.error(f"Error fetching all significant variants: {e}")
            return {"error": str(e)}

def stream_filtered_variants(path, neg_log_cutoff="5e-8"):
    # Ensure cutoff is a string representing a number
    try:
        float_cutoff = float(neg_log_cutoff)
    except ValueError:
        raise ValueError(f"Invalid cutoff value: {neg_log_cutoff}")
    # Use subprocess with argument lists to avoid shell injection
    zcat_proc = subprocess.Popen(['zcat', path], stdout=subprocess.PIPE)
    awk_script = 'NR > 1 && $6 >= cutoff'
    awk_proc = subprocess.Popen(['awk', f'-v', f'cutoff={float_cutoff}', awk_script], stdin=zcat_proc.stdout, stdout=subprocess.PIPE, text=True)
    zcat_proc.stdout.close()  # Allow zcat_proc to receive a SIGPIPE if awk_proc exits.
    for line in awk_proc.stdout:
        yield line.strip().split("\t")
    awk_proc.stdout.close()
    awk_proc.wait()