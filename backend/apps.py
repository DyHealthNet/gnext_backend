from django.apps import AppConfig
from decouple import config
import pandas as pd
from typing import Dict, Tuple

from backend.utils.pheno_cache import get_pheno_df


class ApiConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'backend'

    def ready(selfself):
        try:
            pheno_file = config("PHENO_FILE")
            _ = get_pheno_df(pheno_file)
        except Exception:
            pass




