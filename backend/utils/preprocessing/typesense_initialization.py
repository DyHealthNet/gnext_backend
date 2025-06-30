import typesense
import gzip
import json
import csv
import pandas as pd

from decouple import config

def initialize_typesense(vep_anno_file, pheno_file, batch_size, api_key, typesense_host, typesense_port):

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

    try:
        client.collections["autocomplete"].delete()
    except Exception:
        pass

    client.collections.create(schema_autocomplete)

    # Importing phenotypes
    pheno_dt = pd.read_csv(pheno_file)

    # Importing phenotypes
    for i, r in pheno_dt.iterrows():

        print("Importing phenotype: ", r['phenocode'])

        # Add phenotype to schema
        doc = {
            'type': 'trait',
            'id': str(r['phenocode']),
            'description': r['description'],
            'external_ref': r['external_id'],
            'category': r['category'],
            'filename': r['filename'],
        }
        client.collections["autocomplete"].documents.create(doc)

    # Import variants
    documents = []
    batch_nr = 1
    header = ""
    with gzip.open(vep_anno_file, 'rt') as vcf:
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
                "id": variant_dict["CHROM"] + "_" + variant_dict["POS"] + "_" + variant_dict["REF"] + "/" + variant_dict["ALT"],
                "description": variant_dict["CHROM"] + "_" + variant_dict["POS"] + "_" + variant_dict["REF"] + "/" + variant_dict["ALT"],
                "external_ref": rs_ids,
                "category": "",
                "filename": "",
            }
            documents.append(doc)
            if len(documents) % batch_size == 0:
                payload = '\n'.join(json.dumps(doc) for doc in documents)
                client.collections["autocomplete"].documents.import_(payload, {'action': 'upsert'})
                print("Batch Nr. ", batch_nr, " imported!")
                documents = []
                batch_nr += 1

        # Import remaining documents
        if documents:
            payload = '\n'.join(json.dumps(doc) for doc in documents)
            client.collections["autocomplete"].documents.import_(payload, {'action': 'upsert'})
            print("Batch Nr. ", batch_nr, " imported!")
            documents = []
            batch_nr += 1
