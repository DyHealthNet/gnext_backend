# Script containing functions for extracting data from GWAS summary statistics
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
    for phenotype in phenotypes:

        # Check if the variant is present in the phenotype file
        gwas_res = extract_GWAS_results_for_variant_and_phenotype(phenotype['filename'],
                                                 phenotype['id'],
                                                 phenotype['category'],
                                                 phenotype['description'],
                                                 variant_id)
        if gwas_res is not None:
            phewas_data.append(gwas_res)
    return phewas_data

def extract_variants_for_range(filename, chr, start, end, pval_cutoff=1.0):
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

    try:
        tabix_file = pysam.TabixFile(norm_filepath)
        columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'pvalue', 'beta', 'stderr_beta',
                   'alt_allele_freq']

        for row in tabix_file.fetch(chr, start - 1, end + 1):
            row = row.split("\t")
            row_dict = dict(zip(columns, row))

            # Convert numeric fields
            pvalue = float(row_dict["pvalue"]) if row_dict["pvalue"] != "." else None
            if pvalue is not None and pvalue > pval_cutoff:
                continue  # skip variants above cutoff

            rows.append(row_dict)
            location.append(int(row_dict["pos"]))
            ref.append(row_dict["ref"])
            alt.append(row_dict["alt"])
            external_ids.append(row_dict["rsid"])
            allele_frequencies.append(
                float(row_dict["alt_allele_freq"]) if row_dict["alt_allele_freq"] != "." else None)

        data = {
            "header": columns,
            "rows": rows,
            "location": location,
            "ref": ref,
            "alt": alt,
            "external_ids": external_ids,
            "allele_frequencies": allele_frequencies,
        }

        return data

    except Exception as e:
        print(f"Error fetching data for {chr}:{start}-{end}: {e}")
        return {"error": str(e)}
