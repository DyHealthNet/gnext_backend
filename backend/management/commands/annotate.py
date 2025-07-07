from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
import subprocess
from decouple import config
from django.conf import settings
from django.core.management import CommandError

logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting annotation of GWAS summary statistics files.")
           self.annotate_variants()
           logger.info("Finished annotation of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Annotation with VEP failed: {e}")
           sys.exit(1)

    @staticmethod
    def annotate_variants():
        genome_build = config("GENOME_BUILD")

        GWAS_norm_dir = settings.GWAS_NORM_DIR
        GWAS_vep_dir = settings.GWAS_VEP_DIR
        os.makedirs(GWAS_vep_dir, exist_ok=True)
        GWAS_vcf_file = settings.GWAS_VCF_FILE
        GWAS_annotated_vcf_file = settings.GWAS_ANNO_VCF_FILE
        window_up = config("MAGMA_WINDOW_UP", cast=int)
        window_down = config("MAGMA_WINDOW_DOWN", cast=int)

        # Generate full VCF file TODO: check how to do this in Docker -> because docker pull and docker run executed
        try:
            subprocess.run(["bash", "backend/utils/preprocessing/bash/generate_full_variants_vcf.sh", GWAS_norm_dir, os.path.join(GWAS_vep_dir, GWAS_vcf_file)], check=True)
        except subprocess.CalledProcessError as e:
            raise CommandError(f"Failed to generate full VCF file: {e}")

        # Download VEP docker image and install cache TODO: check how to do this in Docker -> because docker pull and docker run executed
        try:
            subprocess.run(["bash", "backend/utils/preprocessing/bash/setup_vep.sh", GWAS_vep_dir, genome_build], check = True)
        except subprocess.CalledProcessError as e:
            raise CommandError(f"Failed to setup VEP: {e}")

        # Run VEP annotation TODO: needs to be changed when backend in Docker -> because docker run executed
        try:
            subprocess.run(["bash", "../../utils/preprocessing/bash/run_vep.sh", GWAS_vcf_file, GWAS_annotated_vcf_file, GWAS_vep_dir, genome_build, window_up, window_down], check=True)
        except subprocess.CalledProcessError as e:
            raise CommandError(f"Failed to run VEP: {e}")