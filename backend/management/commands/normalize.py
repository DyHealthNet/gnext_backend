from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
from decouple import config
import os
from backend.utils.preprocessing.zorp.zorp import sniffers
from backend.utils.preprocessing.zorp.zorp import parsers
from django.conf import settings
from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from backend.utils.preprocessing.magma.magma import get_bool

logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting normalization of GWAS summary statistics files.")
           self.normalize_GWAS_stats_files()
           logger.info("Finished normalization of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Normalization or LMDB creations of GWAS files failed: {e}")
           sys.exit(1)


    @staticmethod
    def normalize_GWAS_stats_files():
        GWAS_dir = config("GWAS_DIR")
        pheno_file = config("PHENO_FILE")

        GWAS_norm_dir = settings.GWAS_NORM_DIR
        os.makedirs(GWAS_norm_dir, exist_ok = True)

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        parser_options = {
            "chrom_col": int(config("CHR_COLUMN")),
            "pos_col": int(config("POS_COLUMN")),
            "ref_col": int(config("REF_COLUMN")),
            "alt_col": int(config("ALT_COLUMN")),
            "pval_col": int(config("PVAL_COLUMN")),
            "is_neg_log_pvalue": get_bool(config("PVAL_NEGLOG10")),
            'beta_col': int(config("BETA_COLUMN")),
            'stderr_beta_col': int(config("SE_COLUMN")),
            'allele_freq_col': int(config("AF_COLUMN")),
            'rsid': None
        }

        logger.info("Parser options: %s", parser_options)

        # Normalize GWAS files
        with ProcessPoolExecutor(max_workers=int(config("MAX_WORKERS"))) as executor:
            futures = [
                executor.submit(Command.process_and_normalize, r, GWAS_dir, GWAS_norm_dir, parser_options)
                for i, r in pheno_dt.iterrows()
            ]
            for future in as_completed(futures):
                filename = future.result()

    @staticmethod
    def process_and_normalize(r, GWAS_dir, GWAS_norm_dir, parser_options):
        in_filepath = os.path.join(GWAS_dir, r['filename'])
        norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(r['filename']))
        norm_filepath = os.path.join(GWAS_norm_dir, norm_filename)

        if os.path.exists(norm_filepath + ".gz"):
            logger.info("Skipping normalization for %s, file already exists.", norm_filepath)
            return r['filename']

        parser = parsers.GenericGwasLineParser(**parser_options)
        reader = sniffers.guess_gwas_generic(in_filepath, parser=parser, skip_errors=True)
        columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'pvalue', 'beta', 'stderr_beta',
                   'alt_allele_freq']
        reader.write(norm_filepath, make_tabix=True, columns=columns)
        logger.info("COMPLETED: Normalization of GWAS file: %s", r['filename'])
        return r['filename']

