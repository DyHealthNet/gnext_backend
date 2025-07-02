# Script containing functions for extracting data from GWAS summary statistics
from backend.utils.typesense_client import get_all_phenotypes_from_typesense
import pysam
from decouple import config
import gzip
from backend.utils.converters import convert_variant_id

ref_column = config("REF_COLUMN")
alt_column = config("ALT_COLUMN")
pval_column = config("PVAL_COLUMN")
se_column = config("SE_COLUMN")
beta_column = config("BETA_COLUMN")
af_column = config("AF_COLUMN")

def extract_columns_from_file(filepath):
    with gzip.open(filepath, "rt") as f:
        header_line = f.readline().strip()
        columns = header_line.split("\t")
        return columns

def extract_GWAS_results_for_variant_and_phenotype(filename, phenocode, pheno_category, pheno_label, variant_id):
    # Extract GWAS results of variant from phenotype file via tabix
    tabix_file = pysam.TabixFile(filename)
    # Get variant information
    chr, pos, ref, alt = convert_variant_id(variant_id)

    columns = extract_columns_from_file(filename)

    try:
        for row in tabix_file.fetch(chr, pos-1, pos):
            row = row.split("\t")
            row = dict(zip(columns, row))
            if row[ref_column] == ref and row[alt_column] == alt:
                # Format for PheWas plot
                return {"id": phenocode,
                        "trait_group": pheno_category,
                        "trait_label": pheno_label,
                        "log_pvalue": float(row[pval_column]),
                        "se": float(row[se_column]),
                        "beta": float(row[beta_column]),
                        "af": float(row[af_column])
                        }

    except Exception as e:
        print("Error fetching data:" + str(e))
    return None # no match found

def extract_phenotype_results_for_variant(variant_id):
    # Iterate over all phenotypes and extract variant from corresponding file
    phenotypes = get_all_phenotypes_from_typesense() # TODO: check if getting the data from typesense is more efficient or just read the CSV again or storing the CSV in memory
    phewas_data = []
    for phenotype in phenotypes:
        phenotype_file = config("GWAS_DIR") + phenotype['filename']
        # Check if the variant is present in the phenotype file
        gwas_res = extract_GWAS_results_for_variant_and_phenotype(phenotype_file,
                                                 phenotype['id'],
                                                 phenotype['category'],
                                                 phenotype['description'],
                                                 variant_id)
        if gwas_res is not None:
            phewas_data.append(gwas_res)
    return phewas_data
