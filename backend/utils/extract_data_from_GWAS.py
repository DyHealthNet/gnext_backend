# Script containing functions for extracting data from GWAS summary statistics
import contextlib

from backend.utils.typesense_client import get_all_phenotypes_from_typesense
import pysam
from decouple import config
import gzip
from backend.utils.converters import convert_variant_id
from django.conf import settings
import os
import logging
import re
from django.conf import settings
import lmdb
import struct, msgpack

import time


logger = logging.getLogger('backend')

def extract_GWAS_results_for_variant_and_phenotype(filename, phenocode, pheno_category, pheno_label, variant_id):
    # Extract GWAS results of variant from phenotype file via tabix
    norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(filename))
    norm_filepath = os.path.join(settings.GWAS_NORM_DIR, norm_filename + ".gz")

    tabix_file = pysam.TabixFile(norm_filepath)
    # Get variant information
    chr, pos, ref, alt = convert_variant_id(variant_id)

    columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'pvalue', 'beta', 'stderr_beta', 'alt_allele_freq']
    try:
        for row in tabix_file.fetch(chr, pos-1, pos):
            row = row.split("\t")
            row = dict(zip(columns, row))
            if row["ref"] == ref and row["alt"] == alt:
                # Format for PheWas plot
                return {"id": phenocode,
                        "trait_group": pheno_category,
                        "trait_label": pheno_label,
                        "log_pvalue": float(row["pvalue"]) if row["pvalue"] != "." else None,
                        "se": float(row["stderr_beta"]) if row["stderr_beta"] != "." else None,
                        "beta": float(row["beta"]) if row["beta"] != "." else None,
                        "af": float(row["alt_allele_freq"]) if row["alt_allele_freq"] != "." else None,
                        }

    except Exception as e:
        print("Error fetching data:" + str(e))
    return None # no match found

def extract_phenotype_results_for_variant(variant_id):
    # Iterate over all phenotypes and extract variant from corresponding file
    phenotypes = get_all_phenotypes_from_typesense() # TODO: check if getting the data from typesense is more efficient or just read the CSV again or storing the CSV in memory
    phewas_data = []
    i = 0
    for phenotype in phenotypes:
        i+=1

        # Check if the variant is present in the phenotype file
        gwas_res = extract_GWAS_results_for_variant_and_phenotype(phenotype['filename'],
                                                 phenotype['id'],
                                                 phenotype['category'],
                                                 phenotype['description'],
                                                 variant_id)
        if gwas_res is not None:
            phewas_data.append(gwas_res)
        if i % 1000 == 0:
            logger.info(f"Extracted {i} phenotypes for variant {variant_id}")
    return phewas_data

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

    lmdb_path = os.path.join(settings.GWAS_NORM_DIR, f"lmdb_{config('VITE_GENOME_BUILD')}") + "/data.mdb"
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
                        key_bytes = struct.pack('I', pos)
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
