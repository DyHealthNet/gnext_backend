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
from backend.utils.typesense_client import get_all_phenotypes_from_typesense
from backend.utils.lmdb_gene_query import LMDBGeneQuery
from backend.utils.lmdb_gene_magma_query import LMDBGeneMAGMAQuery

logger = logging.getLogger("backend")

# Cache for NF params
_nf_params_cache = None

def get_nf_params():
    """Load NF params file once and cache in memory."""
    global _nf_params_cache
    if _nf_params_cache is None:
        with open(settings.NF_PARAM_FILE, "r") as f:
            _nf_params_cache = json.load(f)
        logger.debug(f"Loaded NF params from {settings.NF_PARAM_FILE}")
    return _nf_params_cache

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
        "neg_log_pvalue": f"{settings.CHR_BGZ_DIR}/chr_{chr}_neg_log_pvalue.tsv.bgz",
        "beta": f"{settings.CHR_BGZ_DIR}/chr_{chr}_beta.tsv.bgz",
        "stderr_beta": f"{settings.CHR_BGZ_DIR}/chr_{chr}_stderr_beta.tsv.bgz",
        "alt_allele_freq": f"{settings.CHR_BGZ_DIR}/chr_{chr}_alt_allele_freq.tsv.bgz",
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
        if metric_data["neg_log_pvalue"][idx] in (".", "NA", "", None):
            logger.info(f"Skipping trait {trait_code} with missing neg_log_pvalue at index {idx}")
            continue

        cat, desc = pheno_map.get(trait_code, ("", trait_code))
        results.append({
            "id": trait_code,
            "x": idx,
            "trait_group": cat,
            "trait_label": desc,
            "neg_log_pvalue": _to_float(metric_data["neg_log_pvalue"][idx]), # for table
            "log_pvalue": _to_float(metric_data["neg_log_pvalue"][idx]), # for plot
            "pvalue": np.power(10, -_to_float(metric_data["neg_log_pvalue"][idx])),
            "beta": _to_float(metric_data["beta"][idx]),
            "stderr_beta": _to_float(metric_data["stderr_beta"][idx]),
            "alt_allele_freq":_to_float(metric_data["alt_allele_freq"][idx])
        })
    return results, min_af, max_af

def extract_variants_for_range(phenocode, chr, start, end, pval_cutoff=1.0, max_rows=10000):
    norm_filename = os.path.join(settings.NORM_DIR, phenocode + ".gz")

    db_handles = {}
    heap = []
    lmdb_gene_file = settings.LMDB_GENE_FILE

    tabix_file = pysam.TabixFile(norm_filename)
    # If nr_samples in norm file -> then nr_samples is another additional column

    columns = ['chrom', 'pos', 'rsid','ref', 'alt', 'neg_log_pvalue',
               'pvalue', 'beta', 'stderr_beta', 'alt_allele_freq']
    # Get NF params from cache
    nf_params = get_nf_params()
    if nf_params.get("n_samples_column", False):
        columns.append('nr_samples')
    col_idx = {name: i for i, name in enumerate(columns)}
    idx_nlog = col_idx['neg_log_pvalue']
    neg_log_cutoff = -math.log10(pval_cutoff)
    try:
        for row in tabix_file.fetch(chr, start - 1, end):
            row = row.split("\t")
            if len(row) != len(columns):
                logger.warning(
                    f"Length of rows ({len(row)}) does not match length of columns ({len(columns)}). Skipping malformed row: {row}")
                continue

            try:
                neg_log_pval = float(row[idx_nlog]) if row[idx_nlog] not in (".", "NA", "") else None
            except (ValueError, TypeError):
                neg_log_pval = None
            if neg_log_pval is None or neg_log_pval < neg_log_cutoff:
                continue

            if len(heap) < max_rows:
                heapq.heappush(heap, (neg_log_pval, row))  # directly store the neg_log_pval
            else:
                if neg_log_pval > heap[0][0]:  # if current value is greater than the smallest in the heap
                    heapq.heapreplace(heap, (neg_log_pval, row))

        rows = [dict(zip(columns, r)) for _, r in heap]
        with LMDBGeneQuery(lmdb_gene_file) as gene_query:
            for r in rows:

                pos = r['pos']
                chr = r['chrom']
                ref = r['ref']
                alt = r['alt']

                r['variant_id'] = f"{chr}_{pos}_{ref}/{alt}"
                r['pvalue'] = float(r['pvalue']) if r['pvalue'] not in (".", "NA", "", None) else None
                r["nearest_genes"] = gene_query.get_genes_for_variant(chr, pos)
                # remove rsid field
                if 'rsid' in r:
                    del r['rsid']
            # remove rsid field from columns
            columns.remove('chrom')
            columns.remove('pos')
            columns.remove('ref')
            columns.remove('alt')
            if 'rsid' in columns:
                columns.remove('rsid')
        data = {
            "header": ['variant_id', 'nearest_genes'] + columns,
            "rows": rows,
        }

        return data

    except Exception as e:
        logger.error(f"Error fetching data for {chr}:{start}-{end}: {e}")
        return {"error": str(e)}

def get_all_sign_variants_cutoff(phenocode, pval_cutoff=5e-8, max_rows=10000):
    norm_filename = os.path.join(settings.NORM_DIR, phenocode + ".gz")

    db_handles = {}
    heap = []
    lmdb_gene_file = settings.LMDB_GENE_FILE

    columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue',
               'pvalue', 'beta', 'stderr_beta', 'alt_allele_freq']
    # Get NF params from cache
    nf_params = get_nf_params()
    if nf_params.get("n_samples_column", False):
        columns.append('nr_samples')
    col_idx = {name: i for i, name in enumerate(columns)}
    idx_nlog = col_idx['neg_log_pvalue']
    neg_log_cutoff = -math.log10(pval_cutoff)

    
    try:
        for row in stream_filtered_variants(norm_filename, neg_log_cutoff=neg_log_cutoff):
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
        with LMDBGeneQuery(lmdb_gene_file) as gene_query:
            for r in rows:
                pos = r['pos']
                chr = r['chrom']
                ref = r['ref']
                alt = r['alt']
                r['variant_id'] = f"{chr}_{pos}_{ref}/{alt}"
                r['pvalue'] = float(r['pvalue']) if r['pvalue'] not in (".", "NA", "", None) else None
                r["nearest_genes"] = gene_query.get_genes_for_variant(chr, pos)
                # remove chr, pos, ref, alt
                del r['chrom']
                del r['pos']
                del r['ref']
                del r['alt']
                # remove rsid field
                if 'rsid' in r:
                    del r['rsid']
            # remove rsid field from columns
            if 'rsid' in columns:
                columns.remove('rsid')
            columns.remove('chrom')
            columns.remove('pos')
            columns.remove('ref')
            columns.remove('alt')
            data = {
                "header": ['variant_id', 'nearest_genes'] + columns,
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

def extract_gene_signals(chr, start, end, strand, gene_id=None, max_rows = 100):
    # Paths to metric files
    metric_files = {
        "neg_log_pvalue": f"{settings.CHR_BGZ_DIR}/chr_{chr}_neg_log_pvalue.tsv.bgz",
        "beta": f"{settings.CHR_BGZ_DIR}/chr_{chr}_beta.tsv.bgz",
        "stderr_beta": f"{settings.CHR_BGZ_DIR}/chr_{chr}_stderr_beta.tsv.bgz",
        "alt_allele_freq": f"{settings.CHR_BGZ_DIR}/chr_{chr}_alt_allele_freq.tsv.bgz",
    }

    # Get NF params from cache
    nf_params = get_nf_params()

    window_up = int(nf_params["window_up"]) * 1000
    window_down = int(nf_params["window_down"]) * 1000

    # Adjust window based on strand
    if strand == "+":
        region_start = max(0, start - window_up)
        region_end = end + window_down
    else:
        region_start = max(0, start - window_down)
        region_end = end + window_up

    # Extract data from all metric files
    data_by_variant = {}

    with pysam.TabixFile(metric_files["neg_log_pvalue"]) as tbx_p:
        # Get trait names from header
        header_line = tbx_p.header[0].lstrip("#")
        trait_names = header_line.split("\t")[5:]

        # Fetch variants
        cand_lines = list(tbx_p.fetch(chr, region_start, region_end))

        if not cand_lines:
            return None

        for line in cand_lines:
            fields = line.split("\t")
            variant_id = fields[4]
            # Convert "." to None, keep the same length as trait_names
            data_by_variant[variant_id] = {
                "chr": fields[0],
                "pos": int(fields[1]),
                "ref": fields[2],
                "alt": fields[3],
                "neg_log_pvalue": [float(x) if x != "." else None for x in fields[5:]]
            }

    # Extract beta values
    with pysam.TabixFile(metric_files["beta"]) as tbx:
        for line in tbx.fetch(chr, region_start, region_end):
            fields = line.split("\t")
            variant_id = fields[4]
            if variant_id in data_by_variant:
                data_by_variant[variant_id]["beta"] = [float(x) if x != "." else None for x in fields[5:]]

    # Extract stderr_beta values
    with pysam.TabixFile(metric_files["stderr_beta"]) as tbx:
        for line in tbx.fetch(chr, region_start, region_end):
            fields = line.split("\t")
            variant_id = fields[4]
            if variant_id in data_by_variant:
                data_by_variant[variant_id]["stderr_beta"] = [float(x) if x != "." else None for x in fields[5:]]

    # Extract alt_allele_freq values
    with pysam.TabixFile(metric_files["alt_allele_freq"]) as tbx:
        for line in tbx.fetch(chr, region_start, region_end):
            fields = line.split("\t")
            variant_id = fields[4]
            if variant_id in data_by_variant:
                data_by_variant[variant_id]["alt_allele_freq"] = [float(x) if x != "." else None for x in fields[5:]]

    # Build long format dataframe
    rows = []
    for variant_id, variant_data in data_by_variant.items():
        # Determine location relative to gene
        variant_pos = variant_data["pos"]
        # Determine location relative to gene (strand-aware)
        if strand == "+":
            # Positive strand: gene goes 5' -> 3' in increasing coordinates
            if variant_pos < start:
                location = "upstream"
                distance = start - variant_pos
            elif variant_pos > end:
                location = "downstream"
                distance = variant_pos - end
            else:
                location = "within_gene"
                distance = 0
        else:
            # Negative strand: gene goes 5' -> 3' in decreasing coordinates
            if variant_pos > end:
                location = "upstream"  # Higher coordinate = upstream on negative strand
                distance = variant_pos - end
            elif variant_pos < start:
                location = "downstream"  # Lower coordinate = downstream on negative strand
                distance = start - variant_pos
            else:
                location = "within_gene"
                distance = 0

        for i, trait_name in enumerate(trait_names):
            rows.append({
                "variant_id": f"{variant_data['chr']}_{variant_data['pos']}_{variant_data['ref']}/{variant_data['alt']}",
                "chr": variant_data["chr"],
                "pos": variant_data["pos"],
                "ref": variant_data["ref"],
                "alt": variant_data["alt"],
                "location": location,
                "distance": distance,
                "trait_name": trait_name,
                "neg_log_pvalue": variant_data["neg_log_pvalue"][i],
                "beta": variant_data.get("beta", [None] * len(trait_names))[i],
                "stderr_beta": variant_data.get("stderr_beta", [None] * len(trait_names))[i],
                "alt_allele_freq": variant_data.get("alt_allele_freq", [None] * len(trait_names))[i],
            })

    signals_df = pd.DataFrame(rows)

    # Filter out rows with None p-values
    signals_df = signals_df[signals_df["neg_log_pvalue"].notna()]

    # Group by trait and get the top variant (highest neg_log_pvalue) for each trait
    signals_df = signals_df.sort_values("neg_log_pvalue", ascending=False)
    signals_df = signals_df.groupby("trait_name", as_index=False).first()

    # Limit by max_rows (top N traits)
    signals_df = signals_df.head(max_rows)

    # Get phenotype descriptions
    phenotypes = get_all_phenotypes_from_typesense()
    phenotypes_df = pd.DataFrame(phenotypes)

    # Merge signals_df and phenotypes_df to get trait descriptions
    signals_df = signals_df.merge(
        phenotypes_df,
        left_on='trait_name',
        right_on='id',
        how='left'
    )

    # Add p-value column
    signals_df['pvalue'] = signals_df['neg_log_pvalue'].apply(lambda x: np.power(10, -x) if pd.notnull(x) else None)

    # Rename signals
    signals_df = signals_df.rename(columns={"variant_id": "top_variant", "trait_name": "trait_id", "label": "trait_label"})
    signals_df = signals_df[["trait_id", "trait_label", "top_variant", "beta", "stderr_beta", "alt_allele_freq","pvalue", "neg_log_pvalue", "location", "distance"]]

    # Rename NaN to ""
    signals_df = signals_df.fillna("")
    
    # Add MAGMA gene p-values if ensg_id is provided
    # Only if MAGMA in .env
    if config("VITE_MAGMA_SHOW", default="false").lower() == "true":
        try:
            magma_bgz_file = settings.GENE_MAGMA_BGZ_FILE
            magma_index_file = settings.GENE_MAGMA_INDEX_FILE

            with LMDBGeneMAGMAQuery(magma_index_file, magma_bgz_file) as query:
                magma_result = query.get_gene_magma_pvalues(gene_id)
            
            if magma_result:
                # Create a dict mapping trait_id -> magma_pvalue
                magma_pval_map = {
                    tp['trait_id']: tp['pvalue'] 
                    for tp in magma_result['trait_pvalues'] 
                    if tp['pvalue'] is not None
                }
                
                # Add MAGMA p-values to signals_df
                signals_df['MAGMA P-value'] = signals_df['trait_id'].map(magma_pval_map)
                
                logger.info(f"Added MAGMA p-values for {len(magma_pval_map)} traits for gene {gene_id}")
            else:
                logger.warning(f"No MAGMA data found for gene {gene_id}")
        
        except Exception as e:
            logger.error(f"Error adding MAGMA p-values for gene {gene_id}: {e}")

    return signals_df

def extract_region_associations(trait_id, chr, start, end):
    """
    Extract association data for LocusZoom by reusing extract_variants_for_range.
    """
    

    # Use existing function with no p-value cutoff
    result = extract_variants_for_range(
        phenocode=trait_id,
        chr=chr,
        start=start,
        end=end,
        pval_cutoff=1.0,  # Get all variants
        max_rows=100000  # Higher limit for visualization
    )
    
    if "error" in result:
        logger.error("Error in extract_variants_for_range: %s", result.get("error"))
        return None
    
    if not result.get("rows"):
        logger.warning("No rows returned for trait_id: %s in region %s:%d-%d", trait_id, chr, region_start, region_end)
        return None

    # Transform to LocusZoom format
    variants = []
    for row in result["rows"]:
        variant = {
            "chromosome": row["chrom"],
            "position": int(row["pos"]),
            "ref_allele": row["ref"],
            "alt_allele": row["alt"],
            "variant": row["chrom"] + ":" + str(row["pos"]) + "_" + row["ref"] + "/" + row["alt"],
            "id": row["chrom"] + "_" + str(row["pos"]) + "_" + row["ref"] + "/" + row["alt"],
            "pvalue": row["pvalue"],
            "log_pvalue": float(row["neg_log_pvalue"]) if row["neg_log_pvalue"] not in (".", None) else None,
            "nearest_genes": row.get("nearest_genes", "")
        }

        # Add optional fields if available
        if row.get("beta") not in (".", None, "NA", ""):
            variant["beta"] = float(row["beta"])
        if row.get("stderr_beta") not in (".", None, "NA", ""):
            variant["se"] = float(row["stderr_beta"])
        if row.get("alt_allele_freq") not in (".", None, "NA", ""):
            variant["alt_allele_freq"] = float(row["alt_allele_freq"])

        variants.append(variant)

    if not variants:
        return None

    logger.info("Extracted %d variants for LocusZoom", len(variants))

    return {
        "data": variants
    }

def extract_gene_magma_pvalues(ensg_id, pval_cutoff=0.05, max_results=100):
    """
    Extract MAGMA gene-based p-values for a specific gene across all traits.
    
    Args:
        ensg_id: Ensembl gene ID (e.g., "ENSG00000123456")
        pval_cutoff: Maximum p-value to include (default: 0.05)
        max_results: Maximum number of results to return (default: 100)
        
    Returns:
        dict with keys:
            - 'ensg_id': Gene ID
            - 'symbol': Gene symbol
            - 'header': List of column names
            - 'rows': List of trait results sorted by p-value
        Returns None if gene not found or error occurs
    """
    try:
        # Get file paths from settings
        magma_bgz_file = settings.GENE_MAGMA_BGZ_FILE
        magma_index_file = settings.GENE_MAGMA_INDEX_FILE
        
        # Check if files exist
        if not os.path.exists(magma_bgz_file):
            logger.error(f"MAGMA BGZ file not found: {magma_bgz_file}")
            return None
        
        if not os.path.exists(magma_index_file):
            logger.error(f"MAGMA index file not found: {magma_index_file}")
            return None
        
        # Query MAGMA p-values
        with LMDBGeneMAGMAQuery(magma_index_file, magma_bgz_file) as query:
            result = query.get_gene_magma_pvalues(ensg_id)
        
        if result is None:
            logger.warning(f"No MAGMA data found for gene {ensg_id}")
            return None
        
        # Get phenotype descriptions from Typesense
        phenotypes = get_all_phenotypes_from_typesense()
        phenotypes_map = {p['id']: p for p in phenotypes}
        
        # Filter and format results
        rows = []
        for trait_pval in result['trait_pvalues']:
            trait_id = trait_pval['trait_id']
            pvalue = trait_pval['pvalue']
            
            # Skip if no p-value or above cutoff
            if pvalue is None or pvalue > pval_cutoff:
                continue
            
            # Get trait info
            trait_info = phenotypes_map.get(trait_id, {})
            trait_label = trait_info.get('label', trait_id)
            trait_group = trait_info.get('category', '')
            
            rows.append({
                'trait_id': trait_id,
                'trait_label': trait_label,
                'trait_group': trait_group,
                'pvalue': pvalue,
                'neg_log_pvalue': -np.log10(pvalue) if pvalue > 0 else None
            })
        
        # Sort by p-value (ascending) and limit results
        rows.sort(key=lambda x: x['pvalue'] if x['pvalue'] is not None else float('inf'))
        rows = rows[:max_results]
        
        logger.info(f"Found {len(rows)} MAGMA associations for gene {ensg_id} (cutoff: {pval_cutoff})")
        
        return {
            'ensg_id': result['ensg_id'],
            'symbol': result['symbol'],
            'header': ['trait_id', 'trait_label', 'trait_group', 'pvalue', 'neg_log_pvalue'],
            'rows': rows
        }
        
    except Exception as e:
        logger.error(f"Error extracting MAGMA p-values for gene {ensg_id}: {e}")
        return None