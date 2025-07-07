from django.core.management.base import BaseCommand
import sys
import traceback
import logging
from decouple import config
import os

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
        try:
            self.check_env_variables()
        except Exception as e:
            traceback.print_exc()
            logger.error(f"Environment file check failed: {e}")
            sys.exit(1)

    @staticmethod
    def check_env_variables():
        """
        Check presence and basic validity of required environment variables using decouple.config.
        """
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
                logger.error(f"{key} is set but path does not exist: {value}")
                sys.exit(1)

            if key in ["CHR_COLUMN", "POS_COLUMN", "REF_COLUMN", "ALT_COLUMN",
                       "PVAL_COLUMN", "SE_COLUMN", "BETA_COLUMN", "AF_COLUMN",
                       "BATCH_SIZE", "MAGMA_WINDOW_UP", "MAGMA_WINDOW_DOWN"]:
                try:
                    int(value)
                except ValueError:
                    logger.error(f"{key} should be an integer but got: {value}")
                    sys.exit(1)

        if missing:
            logger.error(f"Missing environment variables: {', '.join(missing)}")
            sys.exit(1)
        else:
            logger.info("All required environment variables are present.")