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
from concurrent.futures import ProcessPoolExecutor, as_completed


logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           if settings.MAGMA_ENABLED:
               logger.debug("Magma is disabled.")
               return
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

        GWAS_magma_result_dir = settings.GWAS_MAGMA_RESULT_DIR
        os.makedirs(GWAS_magma_result_dir, exist_ok=True)

        GWAS_magma_norm_dir = os.path.join(settings.GWAS_MAGMA_DIR, "input_GWAS_norm")
        gene_annot = os.path.join(settings.GWAS_MAGMA_DIR, settings.GWAS_ANNO_MAGMA_FILE)
        LD_path = config('MAGMA_LD_REF_DIR')
        pheno_file = config('PHENO_FILE')
        pheno_dt = pd.read_csv(pheno_file)

        n_samples = config('N_SAMPLES')
        magma_model = config('MAGMA_MODEL')

        with ProcessPoolExecutor(max_workers=int(config("MAX_WORKERS"))) as executor:
            futures = [
                executor.submit(
                    Command.process_magma_run,
                    r,
                    GWAS_magma_norm_dir,
                    GWAS_magma_result_dir,
                    magma,
                    LD_path,
                    gene_annot,
                    n_samples,
                    magma_model
                )
                for i, r in pheno_dt.iterrows() # 30820 has by far the fastest magma run (under 2min) -> good for testing
            ]
            for future in as_completed(futures):
                filename = future.result()
                logger.info(f"COMPLETED: Generated magma file: %s", filename)

    def process_magma_run(pheno, GWAS_magma_norm_dir, GWAS_magma_result_dir, magma, LD_path, gene_annot, n_samples, magma_model):
        import logging
        logger = logging.getLogger("backend")

        sample_file = os.path.join(GWAS_magma_norm_dir,  pheno['filename'].split(".")[0] + ".txt")
        magma_file = os.path.join(GWAS_magma_result_dir,  pheno['filename'].split(".")[0])

        if not os.path.exists(magma_file + ".genes.out"):

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

        return pheno['filename']