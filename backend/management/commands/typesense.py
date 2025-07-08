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

           # Config variables
           api_key = config('TYPESENSE_KEY')
           typesense_host = config('TYPESENSE_HOST')
           typesense_port = config('TYPESENSE_PORT')
           batch_size = int(config("BATCH_SIZE"))
           pheno_file = config('PHENO_FILE')

           GWAS_annotated_vcf_file = os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)

           # Initialize the typesense client
           logger.info("Starting filling of typesense.")
           self.fill_typesense(GWAS_annotated_vcf_file, pheno_file, batch_size, api_key, typesense_host, typesense_port)
           logger.info("Finished filling of typesense!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Typesense initialization failed: {e}")
           sys.exit(1)

    @staticmethod
    def fill_typesense(GWAS_annotated_vcf_file, pheno_file, batch_size, api_key, typesense_host, typesense_port):
        """
        Fill Typesense with the provided parameters.
        """
        client = typesense.Client({
            'api_key': api_key,
            'nodes': [{
                'host': typesense_host,
                'port': typesense_port,
                'protocol': 'http'
            }],
            'connection_timeout_seconds': 10
        })

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        # Importing phenotypes (here don't check whether phenotype already in typesense, because upsert does this)
        for i, r in pheno_dt.iterrows():
            logger.info(f"Importing phenotype to typesense: {r['phenocode']}")

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
