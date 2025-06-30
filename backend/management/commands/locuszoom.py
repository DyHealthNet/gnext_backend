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
           logger.info("Starting generation of Manhanttan and QQ files.")
           self.preprocess_GWAS_stats_files()
           logger.info("Finished generation of Manhattan and QQ files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Generation of Manhattan and QQ files failed: {e}")
           sys.exit(1)

    @staticmethod
    def preprocess_GWAS_stats_files():
        GWAS_dir = env("GWAS_DIR")
        pheno_file = env("PHENO_FILE")

        GWAS_norm_dir = os.path.join(GWAS_dir, "GWAS_stats_norm")
        os.makedirs(GWAS_norm_dir, exist_ok = True)

        GWAS_manhattan_dir = os.path.join(GWAS_dir, "GWAS_manhattan")
        os.makedirs(GWAS_manhattan_dir, exist_ok=True)

        GWAS_qq_dir = os.path.join(GWAS_dir, "GWAS_qq")
        os.makedirs(GWAS_qq_dir, exist_ok=True)

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        # Normalize GWAS files
        for i, r in pheno_dt.iterrows():
            # Generate Manhattan and QQ JSON file
            manhattan_filepath = os.path.join(GWAS_manhattan_dir, r['file_name'])
            qq_filepath = os.path.join(GWAS_qq_dir, r['file_name'])
            norm_filepath = os.path.join(GWAS_norm_dir, r['file_name'])
            generate_manhattan_qq_json(norm_filepath, manhattan_filepath, qq_filepath)


