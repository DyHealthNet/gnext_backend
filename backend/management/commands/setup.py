from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import os
from django.core import management
import subprocess
from django.core.management.base import CommandError
from django.conf import settings
from decouple import config
import typesense

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting setup of VEP, Typesense, and MAGMA.")
           # Setup VEP -> check if docker pulled and cache installed
           self.setup_VEP()
           # Setup Typesense -> check if typesense running and if schema exists
           self.setup_typesense()
           # Setup MAGMA -> check if executable downloaded
           self.setup_magma()
           logger.info("Finished setup of VEP, Typesense, and MAGMA!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Setup failed: {e}")
           sys.exit(1)

    @staticmethod
    def setup_VEP():
        GWAS_vep_dir = settings.GWAS_VEP_DIR
        genome_build = config("GENOME_BUILD")

        # Check if VEP cache is already downloaded -> if not -> pull docker image of VEP and download cache
        try:
            subprocess.run(["bash", "backend/utils/preprocessing/bash/setup_vep.sh", GWAS_vep_dir, genome_build],
                           check=True)
        except subprocess.CalledProcessError as e:
            raise CommandError(f"Failed to setup VEP: {e}")

    @staticmethod
    def setup_typesense():
        # Check if Typesense is running
        try:
            subprocess.run(["bash", "backend/utils/preprocessing/bash/setup_typesense.sh"],
                           check=True)
        except subprocess.CalledProcessError as e:
            raise CommandError(f"Failed to setup VEP: {e}")

        # Config variables
        api_key = config('TYPESENSE_KEY')
        typesense_host = config('TYPESENSE_HOST')
        typesense_port = config('TYPESENSE_PORT')

        # Check if Autocomplete schema exists -> if not create it
        client = typesense.Client({
            'api_key': api_key,
            'nodes': [{
                'host': typesense_host,
                'port': typesense_port,
                'protocol': 'http'
            }],
            'connection_timeout_seconds': 10
        })

        schema_autocomplete = {
            "name": "autocomplete",
            "fields": [
                {"name": "type", "type": "string", "facet": True},
                {"name": "id", "type": "string"},
                {"name": "description", "type": "string"},
                {"name": "external_ref", "type": "string"},
                {"name": "category", "type": "string"},
                {"name": "filename", "type": "string"},
            ]
        }

        # Check if collection exists, if not create it
        collections = client.collections.retrieve()
        if "autocomplete" not in [collection['name'] for collection in collections]:
            try:
                client.collections.create(schema_autocomplete)
                logger.info("Typesense autocomplete schema created !")
            except Exception as e:
                raise CommandError(f"Failed to create Typesense schema: {e}")
        else:
            # If the collection exists, we can skip this step
            logger.info("Typesense autocomplete schema already exists, skipping creation.")


    @staticmethod
    def setup_magma(): # TODO: Bastienne
        pass