from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import environ
import os
from backend.utils.preprocessing.typesense_initialization import initialize_typesense
from backend.utils.preprocessing.normalize_GWAS_stats_files import normalize_gwas_stats_file
from backend.utils.preprocessing.generate_manhattan_qq_files import generate_manhattan_qq_json
import subprocess

logger = logging.getLogger("backend")

environ.Env.read_env()

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:

           # Config variables
           api_key = config('TYPESENSE_KEY')
           typesense_host = config('TYPESENSE_HOST')
           typesense_port = config('TYPESENSE_PORT')
           batch_size = int(config("BATCH_SIZE"))
           pheno_file = config('PHENO_FILE')

           GWAS_dir = env("GWAS_DIR")
           GWAS_vcf_dir = os.path.join(GWAS_dir, "GWAS_vcf")
           os.makedirs(GWAS_vcf_dir, exist_ok=True)
           GWAS_annotated_vcf_file = os.path.join(GWAS_vcf_dir, "annotated_full_variants.vcf.bgz")

           # Initialize the typesense client
           logger.info("Starting typesense initialization.")
           initialize_typesense(GWAS_annotated_vcf_file, pheno_file, batch_size, api_key, typesense_host, typesense_port)
           logger.info("Finished typesense initialization!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Typesense initialization failed: {e}")
           sys.exit(1)

