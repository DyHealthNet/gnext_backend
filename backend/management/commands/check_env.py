from django.core.management.base import BaseCommand
from backend.utils.preprocessing.magma.magma import get_bool, read_magma_config
import logging
from decouple import config
import os
from django.core.management import CommandError
from decouple import UndefinedValueError
import pandas as pd

from dyhealthnetlight.settings import MAGMA_ENABLED

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
        required_keys = [
            "VITE_TYPESENSE_KEY",
            "VITE_TYPESENSE_HOST",
            "VITE_TYPESENSE_PORT",
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
            "MAGMA_ENABLED", # Can be omitted as it's individually tested in the next step
        ]

        missing = []

        try:
            magma_enabled = get_bool(config("MAGMA_ENABLED"))
            if magma_enabled == None:
                raise CommandError(f"MAGMA_ENABLED must be set to either True/False, 1/0, yes/no or on/off")

            if magma_enabled:
                required_keys += [
                    "MAGMA_LD_REF",
                    "MAGMA_CONFIG_FILE",
                    "N_SAMPLES"
                ]
            else:
                logger.info(f"MAGMA is currently disabled. All MAGMA-related preprocessing steps will be skipped. "
    f"To enable MAGMA analysis for your dataset, set MAGMA_ENABLED to True in your .env file.")
        except UndefinedValueError:
            missing.append("MAGMA_ENABLED")



        file_checks = {"PHENO_FILE":"file", "GWAS_DIR":"dir", "MAGMA_LD_REF":"dir", "MAGMA_CONFIG_FILE":"file"}
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
                       "BATCH_SIZE", "MAGMA_WINDOW_UP", "MAGMA_WINDOW_DOWN",
                       "N_SAMPLES"]:
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

        # Check if output directory set, if so, check if it exists
        output_dir = config("OUTPUT_DIR", default=config("GWAS_DIR"))
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir, exist_ok=True)

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

        if MAGMA_ENABLED:
            mconfig_rows = read_magma_config(config("MAGMA_CONFIG_FILE"))
            logger.info("MAGMA config file is valid.")

