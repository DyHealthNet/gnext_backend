from django.core.management.base import BaseCommand, CommandError
import sys
import pandas as pd
import traceback
import logging
from decouple import config
import os
import typesense
import gzip
import json
from django.conf import settings
import subprocess
import time
import requests

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
        try:
            # Setup Typesense -> check if typesense running and if schema exists
            self.setup_typesense()

            # Get parameters from .env
            batch_size = int(config("BATCH_SIZE"))
            pheno_file = config("PHENO_FILE")
            # Initialize the typesense client
            logger.info("Starting filling of Typesense.")
            self.fill_typesense(settings.ANNO_VCF_FILE, pheno_file, batch_size)
            logger.info("Finished filling of Typesense!")
        except Exception as e:
            traceback.print_exc()
            logger.error(f"Typesense initialization failed: {e}")
            sys.exit(1)

    # -------------------------------------------------------------------------
    # Setup Typesense
    # -------------------------------------------------------------------------
    @staticmethod
    def setup_typesense():
        try:
            CONTAINER_NAME = config("VITE_TYPESENSE_HOST")
            PORT = config("VITE_TYPESENSE_PORT")
            API_KEY = config("VITE_TYPESENSE_KEY")
            VOLUME_DIR = config("TYPESENSE_VOLUME_DIR")

            subprocess.run(
                ["bash", "backend/utils/preprocessing/setup_typesense.sh", CONTAINER_NAME, PORT, API_KEY, VOLUME_DIR],
                check=True
            )

            # Wait for Typesense to become ready
            Command.wait_for_typesense("localhost", PORT, API_KEY)

        except subprocess.CalledProcessError as e:
            raise CommandError(f"Failed to setup Typesense: {e}")

        # Initialize client
        client = typesense.Client(settings.TYPESENSE_CONFIG)

        # Define schema
        schema_autocomplete = {
            "name": "autocomplete",
            "fields": [
                {"name": "id", "type": "string"},
                {"name": "type", "type": "string", "facet": True},
                {"name": "label", "type": "string"},
                {"name": "description", "type": "string"},
                {"name": "external_ref", "type": "string"},
                {"name": "category", "type": "string"},
                {"name": "filename", "type": "string"},
                {"name": "nr_samples", "type": "int32"},
            ],
        }

        # Create collection if missing
        collections = client.collections.retrieve()
        if "autocomplete" not in [collection["name"] for collection in collections]:
            try:
                client.collections.create(schema_autocomplete)
                logger.info("Typesense autocomplete schema created!")
            except Exception as e:
                raise CommandError(f"Failed to create Typesense schema: {e}")
        else:
            logger.info("Typesense autocomplete schema already exists, skipping creation.")

    # -------------------------------------------------------------------------
    # Wait for Typesense health endpoint
    # -------------------------------------------------------------------------
    @staticmethod
    def wait_for_typesense(host: str, port: str, api_key: str, timeout: int = 200):
        url = f"http://{host}:{port}/health"
        headers = {"X-TYPESENSE-API-KEY": api_key}

        logger.info(f"Waiting for Typesense at {url} ...")
        start = time.time()
        while time.time() - start < timeout:
            try:
                r = requests.get(url, headers=headers, timeout=200)
                if r.ok:
                    logger.info("✅ Typesense is ready.")
                    return True
            except requests.exceptions.RequestException:
                pass
            time.sleep(2)
        raise TimeoutError("Typesense did not become ready in time.")

    # -------------------------------------------------------------------------
    # Reset Typesense Collection
    # -------------------------------------------------------------------------
    @staticmethod
    def reset_collection_if_needed(client, collection_name):
        try:
            results = client.collections[collection_name].documents.search({
                "q": "*",
                "query_by": "description",
                "per_page": 1
            })

            if results["found"] > 0:
                client.collections[collection_name].delete()
                schema_autocomplete = {
                    "name": "autocomplete",
                    "fields": [
                        {"name": "id", "type": "string"},
                        {"name": "type", "type": "string", "facet": True},
                        {"name": "label", "type": "string"},
                        {"name": "description", "type": "string"},
                        {"name": "external_ref", "type": "string"},
                        {"name": "category", "type": "string"},
                        {"name": "filename", "type": "string"},
                        {"name": "nr_samples", "type": "int32"},
                    ],
                }
                client.collections.create(schema_autocomplete)
            else:
                logger.info(f"No existing documents found in {collection_name}.")
        except Exception as e:
            logger.error(f"Error resetting collection {collection_name}: {e}")

    # -------------------------------------------------------------------------
    # Fill Typesense
    # -------------------------------------------------------------------------
    @staticmethod
    def fill_typesense(GWAS_annotated_vcf_file, pheno_file, batch_size):
        client = typesense.Client(settings.TYPESENSE_CONFIG)

        Command.reset_collection_if_needed(client, "autocomplete")

        ################### Import phenotypes ###################
        pheno_dt = pd.read_csv(pheno_file)
        for i, r in pheno_dt.iterrows():
            logger.info(f"Importing phenotype to typesense: {r['phenocode']}")
            if "external_id" not in r:
                r["external_id"] = ""
            if "nr_samples" not in r:
                r["nr_samples"] = 0

            doc = {
                "type": "trait",
                "id": str(r["phenocode"]),
                "label": str(r["phenocode"]),
                "description": r["description"],
                "external_ref": r["external_id"],
                "category": r["category"],
                "filename": r["filename"],
                "nr_samples": int(r["nr_samples"]),
            }
            client.collections["autocomplete"].documents.upsert(doc)

        ################### Import genes ###################
        logger.info(f"Reading gene file: {settings.GENE_FILE}")
        mapped_genes = pd.read_csv(settings.GENE_FILE, sep="\t")
        logger.info(f"Loaded {len(mapped_genes)} genes")
        
        documents = []
        for i, r in mapped_genes.iterrows():
            doc = {
                "type": "gene",
                "id": str(r["ensg_id"]),  # Keep ensg_id as unique ID
                "label": str(r["symbol"]),  # Display symbol as label
                "description": f"Chr {r['chr']}: {r['start']}-{r['end']} ({r['strand']})",
                "external_ref": str(r["symbol"]),  # Store symbol in external_ref for routing
                "category": "",
                "filename": "",
                "nr_samples": 0,
            }
            documents.append(doc)
        
        logger.info(f"Prepared {len(documents)} gene documents for import")
        payload = "\n".join(json.dumps(doc) for doc in documents)
        result = client.collections["autocomplete"].documents.import_(payload, {"action": "upsert"})
        logger.info(f"Genes imported to typesense! Result: {result}")

        ################### Import variants ###################
        # First, parse the annotation columns from the VCF header
        anno_columns = []
        with gzip.open(GWAS_annotated_vcf_file, "rt") as vcf:
            for line in vcf:
                if line.startswith("##INFO") and "Format:" in line:
                    anno_columns = line.strip().split("|")
                    anno_columns[0] = anno_columns[0].split("Format: ")[1].split('"')[0]
                    anno_columns[-1] = anno_columns[-1].strip('">')
                    break

        # Now process the variants
        documents = []
        batch_nr = 1
        header = []

        with gzip.open(GWAS_annotated_vcf_file, "rt") as vcf:
            for line in vcf:
                if line.startswith("#CHROM"):
                    header = [h.replace("#", "") for h in line.strip().split("\t")]
                    continue
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                variant_dict = dict(zip(header, fields))

                # Parse INFO field into structured dictionaries
                info = variant_dict["INFO"].replace("CSQ=", "").split(",")
                info_dict_list = [dict(zip(anno_columns, i.split("|"))) for i in info]

                # Extract unique rs_ids from the parsed dictionaries
                rs_ids = list(
                    set([d.get("Existing_variation", "") for d in info_dict_list if d.get("Existing_variation", "")]))
                rs_ids = ", ".join(rs_ids)

                doc = {
                    "type": "variant",
                    "id": f"{variant_dict['CHROM']}_{variant_dict['POS']}_{variant_dict['REF']}/{variant_dict['ALT']}",
                    "label": f"{variant_dict['CHROM']}_{variant_dict['POS']}_{variant_dict['REF']}/{variant_dict['ALT']}",
                    "description": f"{variant_dict['CHROM']}_{variant_dict['POS']}_{variant_dict['REF']}/{variant_dict['ALT']}",
                    "external_ref": rs_ids,
                    "category": "",
                    "filename": "",
                    "nr_samples": 0,
                }
                documents.append(doc)

                if len(documents) % batch_size == 0:
                    payload = "\n".join(json.dumps(doc) for doc in documents)
                    client.collections["autocomplete"].documents.import_(payload, {"action": "upsert"})
                    logger.info(f"Batch Nr. {batch_nr} imported to typesense!")
                    documents.clear()
                    batch_nr += 1

            # Import remaining documents
            if documents:
                payload = "\n".join(json.dumps(doc) for doc in documents)
                client.collections["autocomplete"].documents.import_(payload, {"action": "upsert"})
                logger.info(f"Batch Nr. {batch_nr} imported to typesense!")

