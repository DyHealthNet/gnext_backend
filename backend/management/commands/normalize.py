from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import environ
import os
from backend.utils.preprocessing.normalize_GWAS_stats_files import normalize_gwas_stats_file
from backend.utils.preprocessing.generate_manhattan_qq_files import generate_manhattan_qq_json
import subprocess

logger = logging.getLogger("backend")

environ.Env.read_env()

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting normalization of GWAS summary statistics files.")
           self.normalize_GWAS_stats_files()
           logger.info("Finished normalization of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Normalization of GWAS files failed: {e}")
           sys.exit(1)

    @staticmethod
    def preprocess_GWAS_stats_files():
        GWAS_dir = env("GWAS_DIR")
        pheno_file = env("PHENO_FILE")

        chr_column = env("CHR_COLUMN")
        pos_column = env("POS_COLUMN")
        ref_column = env("REF_COLUMN")
        alt_column = env("ALT_COLUMN")
        pval_column = env("PVAL_COLUMN")
        se_column = env("SE_COLUMN")
        beta_column = env("BETA_COLUMN")
        af_column = env("AF_COLUMN")
        pval_neglog10 = env("PVAL_NEGLOG10")

        genome_build = env("GENOME_BUILD")

        GWAS_norm_dir = os.path.join(GWAS_dir, "GWAS_stats_norm")
        os.makedirs(GWAS_norm_dir, exist_ok = True)

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        # Normalize GWAS files
        for i, r in pheno_dt.iterrows():
            in_filepath = os.path.join(GWAS_dir, r['filename'])
            # Normalize GWAS file
            norm_filepath = os.path.join(GWAS_norm_dir, r['file_name'])
            normalize_gwas_stats_file(in_filepath, norm_filepath)

