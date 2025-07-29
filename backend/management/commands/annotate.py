from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
import subprocess
from decouple import config
from django.conf import settings
from django.core.management import CommandError
from backend.utils.preprocessing.magma.magma import read_magma_config

logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       logger.info("Starting annotation of GWAS summary statistics files.")
       self.annotate_variants()
       logger.info("Finished annotation of GWAS summary statistics files!")

    @staticmethod
    def annotate_variants():
        genome_build = config("GENOME_BUILD")

        GWAS_norm_dir = settings.GWAS_NORM_DIR
        GWAS_vep_dir = settings.GWAS_VEP_DIR
        os.makedirs(GWAS_vep_dir, exist_ok=True)
        GWAS_vcf_file = settings.GWAS_VCF_FILE
        GWAS_annotated_vcf_file = settings.GWAS_ANNO_VCF_FILE
        window_up, window_down = read_magma_config(config("MAGMA_CONFIG_FILE"), return_max_window=True)

        # TODO: currently not generating the VCF file even if new phenotypes have been added

        # Generate full VCF file
        # Check if GWAS_annotated_vcf_file already exists
        if os.path.exists(os.path.join(GWAS_vep_dir, GWAS_vcf_file)):
            logger.info("Skipping VCF generation, because VCF file already exists.")
        else:
            try:
                subprocess.run(["bash", "backend/utils/preprocessing/bash/generate_full_variants_vcf.sh", GWAS_norm_dir, os.path.join(GWAS_vep_dir, GWAS_vcf_file)], check=True)
                logger.info("COMPLETED: Generation of VCF file!")
            except subprocess.CalledProcessError as e:
                raise CommandError(f"Failed to generate full VCF file: {e}")

        # Run VEP annotation
        if os.path.exists(os.path.join(GWAS_vep_dir, GWAS_annotated_vcf_file)):
            logger.info("Skipping VCF annotation, because annotated VCF file already exists.")
        else:
            try:
                subprocess.run(["bash", "backend/utils/preprocessing/bash/run_vep.sh", GWAS_vcf_file, GWAS_annotated_vcf_file, GWAS_vep_dir, genome_build, str(window_up), str(window_down)], check=True)
                logger.info("COMPLETED: Annotation of VCF file!")
            except subprocess.CalledProcessError as e:
                raise CommandError(f"Failed to run VEP: {e}")