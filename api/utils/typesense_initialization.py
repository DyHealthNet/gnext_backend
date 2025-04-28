import typesense
import gzip
import json
import csv

client = typesense.Client({
    'api_key': 'xyz',
    'nodes': [{
        'host': 'localhost',
        'port': '8108',
        'protocol': 'http'
    }],
    'connection_timeout_seconds': 2
})

schema = {
    'name': 'GWAS_stats',
    'fields': [
        {'name': '#CHROM', 'type': 'string'},
        {'name': 'BEG', 'type': 'string'},
        {'name': 'END', 'type': 'string'},
        {'name': 'MARKER_ID', 'type': 'string'},
        {'name': 'MAF', 'type': 'string'},
        {'name': 'PVALUE', 'type': 'float'},
        {'name': 'BETA', 'type': 'string'}
    ],
    'default_sorting_field': 'PVALUE',
}

try:
    client.collections['GWAS_stats'].delete()
except Exception:
    pass

client.collections.create(schema)

documents = []

# get pwd
file_path = "../../data/pheno.3.1.epacts.gz"

def safe_float(val):
    return float(val) if val != 'NA' else None

with gzip.open(file_path, 'rt') as f:
    reader = csv.DictReader(f, delimiter = '\t')
    for row in reader:
        # Only take #CHROM, BEG, END, MARKER_ID, MAF, PVALUE, BETA
        row = {
            '#CHROM': row['#CHROM'],
            'BEG': row['BEG'],
            'END': row['END'],
            'MARKER_ID': row['MARKER_ID'],
            'MAF': row['MAF'],
            'PVALUE': float(row['PVALUE']) if row['PVALUE'] != 'NA' else -9999,
            'BETA': row['BETA']
        }

        documents.append(row)

    # Bulk import
    payload = '\n'.join(json.dumps(doc) for doc in documents)

    response = client.collections['GWAS_stats'].documents.import_(
        payload,
        {'action': 'upsert'}
    )

    stats = client.collections['GWAS_stats'].retrieve()
    print("Document count:", stats['num_documents'])

    response = client.collections['GWAS_stats'].documents.search({
        'q': '1:787399_G/T',
        'query_by': 'MARKER_ID',
        'per_page': 100
    })

    for hit in response['hits']:
        print(hit['document'])