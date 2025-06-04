import setuptools.dist

from api.utils.locuszoom.processors import generate_manhattan, generate_qq
from api.utils.zorp.zorp import parsers, sniffers
import api.utils.locuszoom.processors as processors
from decouple import config
import pandas as pd
import os

# Importing phenotypes
pheno_file = config('PHENO_FILE')
pheno_dt = pd.read_csv(pheno_file)

dir_path = config('GWAS_DIR')

parser_options = {
    "chrom_col": 1,
    "pos_col": 2,
    "ref_col": 3,
    "alt_col": 4,
    "pval_col": 8,
    "is_neg_log_pvalue": True,
    'beta': 6,
    'stderr_beta': 7,
    'alt_allele_freq': 5,
    'rsid': None
}

parser = parsers.GenericGwasLineParser(**parser_options)


# Make directories for normalized GWAS stats, Manhattan plots, and QQ plots
os.makedirs(dir_path + "GWAS_stats_norm/", exist_ok=True)
os.makedirs(dir_path + "GWAS_manhattan/", exist_ok=True)
os.makedirs(dir_path + "GWAS_qq/", exist_ok=True)

# Importing phenotypes
for i, r in pheno_dt.iterrows():
    print("Importing phenotype: ", r['phenocode'])
    GWAS_file = r['filename']
    sample_file = dir_path + GWAS_file
    norm_gwas_file = dir_path + "/GWAS_stats_norm/" + GWAS_file
    manhattan_file = dir_path + "/GWAS_manhattan/" + GWAS_file.replace(".tsv.bgz", "_manhattan.json")
    qq_file = dir_path + "/GWAS_qq/" + GWAS_file.replace(".tsv.bgz", "_qq.json")

    reader = sniffers.guess_gwas_generic(sample_file, parser = parser, skip_errors = True)

    status = processors.normalize_contents(
              reader,
               norm_gwas_file,
               'GRCh37',
               debug_mode=True
    )
    print("Normalization successful!")

    generate_manhattan(build = "GRCh37", in_filename = norm_gwas_file + ".gz", out_filename = manhattan_file)
    print("Manhattan JSON generation successful!")

    generate_qq(in_filename = norm_gwas_file + ".gz", out_filename = qq_file)
    print("QQ JSON generation successful!")

