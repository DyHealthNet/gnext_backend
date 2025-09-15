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
from backend.utils.pheno_cache import get_pheno_df, get_pheno_map

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