from django.core.management.base import BaseCommand
import logging
from decouple import config
import os
from django.core.management import CommandError
from decouple import UndefinedValueError
import pandas as pd

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
        required_keys = [
            "TYPESENSE_KEY",
            "TYPESENSE_HOST",
            "TYPESENSE_PORT",
            "VITE_API_URL",
            "VITE_STUDY_NAME",
            "VITE_GENOME_BUILD",
            "FRONTEND_URL",
            "PHENO_FILE",
            "GWAS_DIR",
            "GENOME_BUILD",
            "CHR_COLUMN",
            "POS_COLUMN",
            "REF_COLUMN",
            "ALT_COLUMN",
            "PVAL_COLUMN",
            "SE_COLUMN",
            "BETA_COLUMN",
            "AF_COLUMN",
            "PVAL_NEGLOG10",
            "BATCH_SIZE",
            "MAGMA_WINDOW_UP",
            "MAGMA_WINDOW_DOWN",
            "MAGMA_LD_REF",
        ]

        missing = []
        file_checks = {"PHENO_FILE":"file", "GWAS_DIR":"dir", "MAGMA_LD_REF":"dir"}
        for key in required_keys:
            try:
                value = config(key)
            except UndefinedValueError:
                missing.append(key)
                continue

            if key in file_checks:
                if file_checks[key] == "file" and not os.path.isfile(value):
                    raise CommandError(f"{key} is set but the file does not exist: {value}")
                if file_checks[key] == "directory" and not os.path.isdir(value):
                    raise CommandError(f"{key} is set but the directory does not exist: {value}")

            if key in ["CHR_COLUMN", "POS_COLUMN", "REF_COLUMN", "ALT_COLUMN",
                       "PVAL_COLUMN", "SE_COLUMN", "BETA_COLUMN", "AF_COLUMN",
                       "BATCH_SIZE", "MAGMA_WINDOW_UP", "MAGMA_WINDOW_DOWN"]:
                try:
                    int(value)
                except ValueError:
                    raise CommandError(f"{key} should be an integer but got: {value}")

            if key == "GENOME_BUILD":
                if value not in ["GRCh37", "GRCh38"]:
                    raise CommandError(f"Invalid genome build: {value}. Expected 'GRCh37' or 'GRCh38'.")

        if missing:
            raise CommandError(f"Missing environment variables: {', '.join(missing)}")
        else:
            logger.info("All required environment variables are present.")


        # Check if all required columns in pheno_file are present
        # TODO: Lisi

        # Check if all GWAS stats files are in the GWAS_DIR
        GWAS_dir = config("GWAS_DIR")
        pheno_file = config("PHENO_FILE")

        # Importing phenotype table
        pheno_dt = pd.read_csv(pheno_file)

        for i, r in pheno_dt.iterrows():
            in_filepath = os.path.join(GWAS_dir, r['filename'])
            if not os.path.isfile(in_filepath):
                raise CommandError(f"GWAS file does not exist: {in_filepath}")

        logger.info("All GWAS summary statistics files are present.")