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
            "VITE_TYPESENSE_KEY",
            "VITE_TYPESENSE_HOST",
            "VITE_TYPESENSE_PORT",
            "VITE_GENOME_BUILD",
            "VITE_HG_BUILD_NUMBER",
            "VITE_STUDY_NAME",
            "PHENO_FILE",
            "NF_DATA_DIR",
            "BATCH_SIZE",
            "VITE_API_URL",
        ]

        missing = []

        file_checks = {"PHENO_FILE":"file", "NF_DATA_DIR":"dir"}
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

            try:
                int(config("BATCH_SIZE"))
            except ValueError:
                raise CommandError(f"BATCH_SIZE should be an integer but got: {config('BATCH_SIZE')}")

            if key == "GENOME_BUILD":
                if value not in ["GRCh37", "GRCh38"]:
                    raise CommandError(f"Invalid genome build: {value}. Expected 'GRCh37' or 'GRCh38'.")

        if missing:
            raise CommandError(f"Missing environment variables: {', '.join(missing)}")
        else:
            logger.info("All required environment variables are present.")


        # Check if all required columns in pheno_file are present
        pheno_dt = pd.read_csv(config("PHENO_FILE"))
        required_pheno_cols = ["phenocode", "description", "nr_samples", "category", "filename"]
        for col in required_pheno_cols:
            if col not in pheno_dt.columns:
                raise CommandError(f"Missing required column '{col}' in phenotype file.")
        logger.info("All required columns are present in the phenotype file.")

        # Check if nextflow directory exists with required files
        nf_data_dir = config("NF_DATA_DIR")
        if not os.path.isdir(nf_data_dir):
            raise CommandError(f"NF_DATA_DIR does not exist or is not a directory: {nf_data_dir}")


