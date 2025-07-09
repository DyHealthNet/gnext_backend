from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
from django.core import management
import subprocess
from django.core.management.base import CommandError
from decouple import config
import typesense
from django.conf import settings
import pandas as pd

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting healthcheck of preprocessed data.")
           # Check if environment variables are still correctly set
           management.call_command("check_env")

           # Check if all files exist
           self.check_files()

           # Check if typesense docker is running and has entries in the autocomplete collection
           container_name = "dyhealthnetlight-typesense"
           self.check_typesense(container_name)

           logger.info("Finished healthcheck of preprocessed data!")


       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Healthcheck failed: {e}")
           sys.exit(1)

    @staticmethod
    def check_files():
        GWAS_dir = config("GWAS_DIR")
        pheno_file = config("PHENO_FILE")

        # Importing phenotype table
        pheno_dt = pd.read_csv(pheno_file)

        def check_file_exists(path, description):
            if not os.path.isfile(path):
                raise CommandError(f"{description} does not exist: {path}")

        for i, r in pheno_dt.iterrows():
            in_filepath = os.path.join(GWAS_dir, r['filename'])
            norm_filepath = os.path.join(settings.GWAS_NORM_DIR, r['filename'].split(".")[0] + ".gz")
            manhattan_filepath = os.path.join(settings.GWAS_MANHATTAN_DIR, r['filename'].split(".")[0] + "_manhattan.json")
            qq_filepath = os.path.join(settings.GWAS_QQ_DIR, r['filename'].split(".")[0] + "_qq.json")
            magma_input_filepath = "" #TODO:Bastienne
            magma_output_filepath = "" #TODO: Bastienne

            check_file_exists(in_filepath, f"Input file for phenotype {r['phenocode']}")
            check_file_exists(norm_filepath, f"Normalized file for phenotype {r['phenocode']}")
            check_file_exists(manhattan_filepath, f"Manhattan plot file for phenotype {r['phenocode']}")
            check_file_exists(qq_filepath, f"QQ plot file for phenotype {r['phenocode']}")
            #check_file_exists(magma_input_filepath, f"MAGMA input file for phenotype {r['phenocode']}") #TODO: Bastienne -> uncomment
            #check_file_exists(magma_output_filepath, f"MAGMA output file for phenotype {r['phenocode']}") #TODO: Bastienne -> uncomment


        logger.info("All phenotype-specific files are present.")

        # Check if annotation files exist
        if not os.path.isfile(os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_VCF_FILE)):
            raise CommandError(f"VCF file for annotation does not exist: {settings.GWAS_VCF_FILE}")
        if not os.path.isfile(os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)):
            raise CommandError(f"VEP annotation file does not exist: {settings.GWAS_ANNO_VCF_FILE}")

        # Check if MAGMA mapping file exists
        #if not os.path.isfile(os.path.join(settings.GWAS_MAGMA_DIR, settings.GWAS_ANNO_MAGMA_FILE)):
        #    raise CommandError(f"MAGMA mapping file does not exist: {settings.GWAS_ANNO_MAGMA_FILE}")


    @staticmethod
    def check_typesense(container_name):

        # Build the command string
        command = f'docker ps --filter "name={container_name}" --filter "status=running" --format "{{{{.Names}}}}"'

        try:
            # Run the command and capture output
            result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True)

            # Get the container names from stdout
            output = result.stdout.strip()

            if container_name in output:
                logger.info("Typesense-server container is running.")
            else:
                raise CommandError("Typesense-server container is NOT running.")

        except subprocess.CalledProcessError as e:
            raise CommandError(f"Error checking for typesense docker container: {e}")

        # Check if typesense has entries in the autocomplete collection

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


        # Check if collection exists, if not create it
        collections = client.collections.retrieve()
        if "autocomplete" not in [collection['name'] for collection in collections]:
            raise CommandError("Typesense autocomplete schema does not exist!")

        # Check if there are entries in the autocomplete collection
        try:
            # We run a search with empty query to get all documents, but per_page=0 to avoid actually returning any
            res = client.collections["autocomplete"].documents.search({
                "q": "*",
                "query_by": "description",  # IMPORTANT: replace with the actual field you use for querying
                "per_page": 0
            })

            count = res["found"]
            if count == 0:
                raise CommandError("Typesense autocomplete collection is empty!")
            else:
                logger.info(f"Typesense autocomplete collection has {count} entries.")

        except Exception as e:
            raise CommandError(f"Error checking typesense autocomplete collection: {e}")
