from django.core.management.base import BaseCommand
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

logger = logging.getLogger("backend")


class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           batch_size = int(config("BATCH_SIZE"))
           pheno_file = config('PHENO_FILE')

           GWAS_annotated_vcf_file = os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)

           # Initialize the typesense client
           logger.info("Starting filling of typesense.")
           self.fill_typesense(GWAS_annotated_vcf_file, pheno_file, batch_size)
           logger.info("Finished filling of typesense!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Typesense initialization failed: {e}")
           sys.exit(1)

    @staticmethod
    def reset_collection_if_needed(client, collection_name):
        try:
            # Check if collection has any documents
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
                        {"name": "type", "type": "string", "facet": True},
                        {"name": "id", "type": "string"},
                        {"name": "description", "type": "string"},
                        {"name": "external_ref", "type": "string"},
                        {"name": "category", "type": "string"},
                        {"name": "filename", "type": "string"},
                    ]
                }
                client.collections.create(schema_autocomplete)
            else:
                logger.info(f"No existing documents found in {collection_name}.")

        except Exception as e:
            logger.error(f"Error resetting collection {collection_name}: {e}")

    @staticmethod
    def fill_typesense(GWAS_annotated_vcf_file, pheno_file, batch_size):
        """
        Fill Typesense with the provided parameters.
        """
        client = typesense.Client(config.TYPESENSE_CONFIG)

        # Remove all entries from the collection
        Command.reset_collection_if_needed(client, "autocomplete")

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        # Importing phenotypes (here don't check whether phenotype already in typesense, because upsert does this)
        for i, r in pheno_dt.iterrows():
            logger.info(f"Importing phenotype to typesense: {r['phenocode']}")

            # Check if external_id exists, if not -> add ""
            if 'external_id' not in r:
                r['external_id'] = ""

            # Add phenotype to schema
            doc = {
                'type': 'trait',
                'id': str(r['phenocode']),
                'description': r['description'],
                'external_ref': r['external_id'],
                'category': r['category'],
                'filename': r['filename'],
            }
            client.collections["autocomplete"].documents.upsert(doc)

        # Import variants
        documents = []
        batch_nr = 1
        header = ""
        with gzip.open(GWAS_annotated_vcf_file, 'rt') as vcf:
            for line in vcf:
                if line.startswith("#CHROM"):
                    header = line.strip().split("\t")
                    header = [h.replace("#", "") for h in header]
                    continue
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                variant_dict = dict(zip(header, fields))
                info = variant_dict["INFO"].split(",")
                rs_ids = list(set([i.split("|")[17] for i in info]))
                rs_ids = ", ".join(rs_ids)
                doc = {
                    "type": "variant",
                    "id": variant_dict["CHROM"] + "_" + variant_dict["POS"] + "_" + variant_dict["REF"] + "/" +
                          variant_dict["ALT"],
                    "description": variant_dict["CHROM"] + "_" + variant_dict["POS"] + "_" + variant_dict["REF"] + "/" +
                                   variant_dict["ALT"],
                    "external_ref": rs_ids,
                    "category": "",
                    "filename": "",
                }
                documents.append(doc)
                if len(documents) % batch_size == 0:
                    payload = '\n'.join(json.dumps(doc) for doc in documents)
                    client.collections["autocomplete"].documents.import_(payload, {'action': 'upsert'})
                    logger.info(f"Batch Nr. {batch_nr} imported to typesense!")
                    documents = []
                    batch_nr += 1

            # Import remaining documents
            if documents:
                payload = '\n'.join(json.dumps(doc) for doc in documents)
                client.collections["autocomplete"].documents.import_(payload, {'action': 'upsert'})
                logger.info(f"Batch Nr. {batch_nr} imported to typesense!")
                documents = []
                batch_nr += 1
