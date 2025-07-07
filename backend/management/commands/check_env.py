from django.core.management.base import BaseCommand
import sys
import traceback
import logging
from decouple import config
import os
from django.core.management import CommandError

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
            "MAGMA_EXEC",
            "MAGMA_LD_REF",
        ]

        missing = []
        file_checks = ["PHENO_FILE", "GWAS_DIR", "MAGMA_LD_REF", "MAGMA_EXEC"]
        for key in required_keys:
            try:
                value = config(key)
            except Exception:
                missing.append(key)
                continue

            if key in file_checks and not os.path.exists(value):
                raise CommandError(f"{key} is set but path does not exist: {value}")

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