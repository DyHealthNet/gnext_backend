import typesense
import gzip
import json
import csv

file_path = "../../data/GWAS_Catalog/30804560-GCST007429-EFO_0004312.h.tsv.gz"

client = typesense.Client({
    'api_key': 'xyz',
    'nodes': [{
        'host': 'localhost',
        'port': '8108',
        'protocol': 'http'
    }],
    'connection_timeout_seconds': 10
})

collection_name_gwas = "gwas_stats"
collection_name_nodes = "nodes"

schema = {
  "name": collection_name_gwas,
  "fields": [
    { "name": "label", "type": "string" },
    { "name": "rsid", "type": "string" },
    { "name": "chr", "type": "string" },
    { "name": "pos", "type": "string" },
    {"name": "effect_allele", "type": "string" },
    {"name": "other_allele", "type": "string" },
    {"name": "p_value", "type": "string" },
    {"name": "beta", "type": "string" },
    {"name": "standard_error", "type": "string" },
  ]
}

try:
    client.collections[collection_name_gwas].delete()

except Exception:
    pass

client.collections.create(schema)


# Stream + buffer from .bgz
BATCH_SIZE = 100000
batch_nr = 1

#batch_gwas = []
batch_nodes = []

with gzip.open(file_path, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        # Add unique ID
        # Only take chr, pos, ref, alt, beta_EUR, se_EUR, neglog10pval_EUR
        # and convert to appropriate types
        #gwas_doc = {
        #    'label': "chr" + row['hm_chrom'] + ":" + row['hm_pos'] + " " + row['hm_effect_allele'] + "|" + row['hm_other_allele'],
        #    'rsid': row['hm_rsid'],
        #    'chr': row['hm_chrom'],
        #    'pos': row['hm_pos'],
        #    'effect_allele': row['hm_effect_allele'],
       #     'other_allele': row['hm_other_allele'],
       #     'p_value': row['p_value'],
       #     'beta': row['beta'],
       #     'standard_error': row['standard_error'],
        #}
        #batch_gwas.append(gwas_doc)

        node_doc = {
            'internal_id': "chr" + row['hm_chrom'] + ":" + row['hm_pos'] + " " + row['hm_effect_allele'] + "|" + row['hm_other_allele'],
            'description': "NA",
            'label': row['hm_rsid'],
            'type': 'variant',
            'gene_symbol': "NA",
        }
        batch_nodes.append(node_doc)

        if len(batch_nodes) >= BATCH_SIZE:
            #docs = '\n'.join(json.dumps(doc) for doc in batch_gwas)
            #client.collections[collection_name_gwas].documents.import_(docs, {"action": "upsert"})

            docs = '\n'.join(json.dumps(doc) for doc in batch_nodes)
            client.collections[collection_name_nodes].documents.import_(docs, {"action": "upsert"})

            print("Batch ", batch_nr, " imported!")
            batch_nodes = []
            batch_nr += 1

# Final flush
if batch_nodes:
    #docs = '\n'.join(json.dumps(doc) for doc in batch_gwas)
    #client.collections[collection_name_gwas].documents.import_(docs, {"action": "upsert"})

    docs = '\n'.join(json.dumps(doc) for doc in batch_nodes)
    client.collections[collection_name_nodes].documents.import_(docs, {"action": "upsert"})

    print("Batch ", batch_nr, " imported!")

print("✅ Import complete.")


