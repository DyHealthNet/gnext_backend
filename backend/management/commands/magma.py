from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
import time
import subprocess
from decouple import config
import pandas as pd
from django.conf import settings

logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting MAGMA execution of GWAS summary statistics files.")
           self.run_MAGMA()
           logger.info("Finished MAGMA execution of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"MAGMA run with VEP failed: {e}")
           sys.exit(1)

    def run_MAGMA(self):
        magma = config('MAGMA_EXEC')
        if not os.path.isfile(magma) or not os.access(magma, os.X_OK):
            download_dir = os.path.join(settings.GWAS_MAGMA_DIR, "magma")
            if not os.path.isfile(os.path.join(download_dir, "magma")) or not os.access(os.path.join(download_dir, "magma"),
                                                                                 os.X_OK):
                logger.error("MAGMA executable not found, please run setup again")
                return

        GWAS_magma_norm_dir = os.path.join(settings.GWAS_MAGMA_DIR, "input_GWAS_norm")
        gene_annot = os.path.join(settings.GWAS_MAGMA_DIR, settings.GWAS_ANNO_MAGMA_FILE)
        LD_path = config('MAGMA_LD_REF_DIR')
        pheno_file = config('PHENO_FILE')
        pheno_dt = pd.read_csv(pheno_file)

        n_samples = config('N_SAMPLES')
        magma_model = config('MAGMA_MODEL')

        for i, r in pheno_dt.iterrows():
            gwas_file = r['filename']
            sample_file = os.path.join(GWAS_magma_norm_dir,  r['filename'].split(".")[0] + ".txt")
            magma_file = os.path.join(settings.GWAS_MAGMA_RESULT_DIR,  r['filename'].split(".")[0])

            if not os.path.exists(magma_file):

                gene_cmd = [
                    magma,
                    '--bfile', LD_path,
                    'synonym-dup=drop-dup',  # optional or change!!
                    '--gene-annot', gene_annot,
                    '--pval', sample_file, f'N={n_samples}',
                    '--gene-model', f'{magma_model}',  # optional or change but this ones best and fastest for summary stats
                    '--genes-only',  # since we do not perform pathway analysis with magma we don't need this
                    '--out', magma_file
                ]

                logger.info(f"[INFO] Running gene analysis for {sample_file}")
                start_time = time.time()
                result = subprocess.run(gene_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                logger.debug("MAGMA stdout:\n%s", result.stdout)
                logger.debug("MAGMA stderr:\n%s", result.stderr)
                if result.returncode != 0:
                    logger.error("MAGMA failed with exit code %d", result.returncode)
                    raise subprocess.CalledProcessError(result.returncode, gene_cmd, output=result.stdout,
                                                        stderr=result.stderr)
                end_time = time.time()
                elapsed_time = end_time - start_time
                logger.debug(f"[DONE] Finished MAGMA run {magma_file} in {elapsed_time:.2f} seconds")