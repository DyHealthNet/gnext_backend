import urllib
import zipfile

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
           logger.info("Starting setup of VEP%s." % (", Typesense and MAGMA" if settings.MAGMA_ENABLED else " and Typesense"))
           # Setup VEP -> check if docker pulled and cache installed
           self.setup_VEP()
           # Setup Typesense -> check if typesense running and if schema exists
           self.setup_typesense()
           # Setup MAGMA -> check if executable downloaded if Magma is enabled
           if settings.MAGMA_ENABLED:
               self.setup_magma()
           logger.info("Finished setup of VEP%s." % (", Typesense and MAGMA" if settings.MAGMA_ENABLED else " and Typesense"))
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
    def setup_magma():
        magma_path = config('MAGMA_EXEC')

        # Check if the executable exists under config path and is runnable
        if os.path.isfile(magma_path) and os.access(magma_path, os.X_OK):
            logger.info("Magma executable exists and is runnable: %s", magma_path)
            return

        # Check if the executable exists under default path and is runnable
        download_dir = os.path.join(settings.GWAS_MAGMA_DIR, "magma")
        os.makedirs(download_dir, exist_ok=True)
        if os.path.isfile(os.path.join(download_dir, "magma")) and os.access(os.path.join(download_dir, "magma"),
                                                                             os.X_OK):
            logger.info("Magma executable exists and is runnable: %s", os.path.join(download_dir, "magma"))
            return

        logger.info("MAGMA executable not found, downloading MAGMA...")
        os.makedirs(download_dir, exist_ok=True)

        zip_path = os.path.join(download_dir, "magma.zip")
        url = "https://vu.data.surfsara.nl/index.php/s/lxDgt2dNdNr6DYt/download"

        try:
            urllib.request.urlretrieve(url, zip_path)
            logger.debug(f"Downloaded MAGMA zip to: {zip_path}")
        except Exception as e:
            logger.error(f"Failed to download MAGMA: {e}")
            return

        try:
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(download_dir)
            logger.debug(f"Extracted MAGMA to: {download_dir}")
        except zipfile.BadZipFile as e:
            logger.error(f"Failed to extract MAGMA zip: {e}")
            return

        # Get the executable
        magma_exe = None
        for root, _, files in os.walk(download_dir):
            for file in files:
                if file == "magma":
                    path = os.path.join(root, file)
                    os.chmod(path, 0o755)  # Ensure it's executable
                    magma_exe = path
                    break
            if magma_exe:
                break

        if magma_exe:
            logger.info(f"Magma executable ready at: {magma_exe}")
            return
        else:
            logger.error("Magma executable not found after extraction.")
            return