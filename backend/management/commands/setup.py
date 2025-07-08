from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import os
from django.core import management

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting setup of VEP, Typesense, and MAGMA.")
           # Setup VEP -> check if docker pulled and cache installed -> TODO: Lisi


           # Setup Typesense -> check if typesense running and if schema exists -> TODO: Lisi


           # Setup MAGMA -> check if executable downloaded -> TODO: Bastienne

           logger.info("Finished setup of VEP, Typesense, and MAGMA!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Setup failed: {e}")
           sys.exit(1)
