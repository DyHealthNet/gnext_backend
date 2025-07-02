from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
import subprocess
from decouple import config

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
        GWAS_dir = config("GWAS_DIR")
        genome_build = config("GENOME_BUILD")

        GWAS_norm_dir = os.path.join(GWAS_dir, "GWAS_stats_norm")

        GWAS_vep_dir = os.path.join(GWAS_dir, "GWAS_vep")
        os.makedirs(GWAS_vep_dir, exist_ok=True)
        GWAS_vcf_file = os.path.join(GWAS_vep_dir, "full_variants.vcf")

        # Generate full VCF file TODO: check how to do this in Docker -> because docker pull and docker run executed
        subprocess.run(["bash", "backend/utils/preprocessing/bash/generate_full_variants_vcf.sh", GWAS_norm_dir, GWAS_vcf_file])

        # Download VEP docker image and install cache TODO: check how to do this in Docker -> because docker pull and docker run executed
        subprocess.run(["bash", "backend/utils/preprocessing/bash/setup_vep.sh", GWAS_vep_dir, genome_build])

        # Run VEP annotation TODO: needs to be changed when backend in Docker -> because docker run executed
        subprocess.run(["bash", "../../utils/preprocessing/bash/run_vep.sh", "full_variants.vcf", "annotated_full_variants.vcf.bgz", GWAS_vep_dir,genome_build])
