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
    'connection_timeout_seconds': 10
})

collection_name = "nodes"

# schema = {
#   "name": collection_name,
#   "fields": [
#     { "name": "label", "type": "string" },
#     { "name": "type", "type": "string" },
#     { "name": "description", "type": "string" },
#     { "name": "internal_id", "type": "string" },
#     {"name": "gene_symbol", "type": "string" }
#   ]
# }
#
# try:
#     client.collections[collection_name].delete()
#
# except Exception:
#     pass
#
# client.collections.create(schema)
#
#
# BATCH_SIZE = 100
# file_path = "../../data/CHRIS_somalogic_descriptive_statistic.tsv.gz"
# documents = []
#
# with gzip.open(file_path, 'rt') as f:
#     reader = csv.DictReader(f, delimiter='\t')
#     for i, row in enumerate(reader, start=1):
#         doc = {
#             'internal_id': row['protein_id'],
#             'description': row['TargetFullName'],
#             'label': row['UniProt'],
#             'type': 'protein',
#             'gene_symbol': row['EntrezGeneSymbol'],
#         }
#         documents.append(doc)
#         if i % BATCH_SIZE == 0:
#             payload = '\n'.join(json.dumps(doc) for doc in documents)
#             client.collections[collection_name].documents.import_(payload, {'action': 'upsert'})
#             print("Batch ", i, " imported!")
#             documents = []
#
#     # Final batch
#     if documents:
#         payload = '\n'.join(json.dumps(doc) for doc in documents)
#         client.collections[collection_name].documents.import_(payload, {'action': 'upsert'})
#
#
#
# print("✅ Import complete.")

print(client.collections[collection_name].retrieve()["num_documents"])
