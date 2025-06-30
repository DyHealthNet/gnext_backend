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
           logger.info("Starting annotation of GWAS summary statistics files.")
           self.preprocess_GWAS_stats_files()
           logger.info("Finished annotation of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Annotation with VEP failed: {e}")
           sys.exit(1)

    @staticmethod
    def preprocess_GWAS_stats_files():
        GWAS_dir = env("GWAS_DIR")

        genome_build = env("GENOME_BUILD")

        GWAS_norm_dir = os.path.join(GWAS_dir, "GWAS_stats_norm")
        os.makedirs(GWAS_norm_dir, exist_ok = True)

        GWAS_vcf_dir = os.path.join(GWAS_dir, "GWAS_vcf")
        os.makedirs(GWAS_vcf_dir, exist_ok=True)
        GWAS_vcf_file = os.path.join(GWAS_vcf_dir, "full_variants.vcf")
        GWAS_annotated_vcf_file = os.path.join(GWAS_vcf_dir, "annotated_full_variants.vcf.bgz")

        # Generate full VCF file
        subprocess.run(["bash", "../../utils/preprocessing/bash/generate_full_variants_vcf.sh", GWAS_norm_dir, GWAS_vcf_file])

        # Run VEP annotation TODO: check VEP cache !!!
        subprocess.run(["bash", "../../utils/preprocessing/bash/run_vep.sh", GWAS_vcf_file, GWAS_annotated_vcf_file, genome_build])
