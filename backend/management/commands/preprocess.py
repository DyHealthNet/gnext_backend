from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import environ
import os
from django.core import management

logger = logging.getLogger("backend")

environ.Env.read_env()

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting preprocessing of GWAS summary statistics files.")
           # Call each management command
           management.call_command("normalize")
           management.call_command("locuszoom")
           management.call_command("annotate")
           management.call_command("typesense")
           logger.info("Finished preprocessing of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Preprocessing failed: {e}")
           sys.exit(1)
