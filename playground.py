import time

import pandas as pd
import os
import re
from decouple import config

from backend.utils.extract_data_from_GWAS import extract_variant_metrics
from dyhealthnetlight.settings import GWAS_QQ_DIR

# filepath = "/storage03/larend/pgwas_data/pgwas_full_phenotype_file.csv"
#
# df = pd.read_csv(filepath)
#
# print(df.head())
#
#
#
# for i, row in df.iterrows():
#     # Check if manhattan and QQ file exist
#     GWAS_DIR = config("GWAS_DIR")
#     OUTPUT_DIR = config("OUTPUT_DIR", default=GWAS_DIR)
#     GWAS_MANHATTAN_DIR = os.path.join(OUTPUT_DIR, "GWAS_manhattan")
#     GWAS_QQ_DIR = os.path.join(OUTPUT_DIR, "GWAS_qq")
#     norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(row['filename']))
#     manhattan_filepath = os.path.join(GWAS_MANHATTAN_DIR, norm_filename + "_manhattan.json")
#     qq_filepath = os.path.join(GWAS_QQ_DIR, norm_filename + "_qq.json")
#     if not os.path.exists(manhattan_filepath):
#         print(f"Manhattan file does not exist: {manhattan_filepath}")
#     if not os.path.exists(qq_filepath):
#         print(f"QQ file does not exist: {qq_filepath}")
#
#

# from django.conf import settings
# import typesense
#
# TYPESENSE_CONFIG = {
#     "nodes": [{
#         "host": config("VITE_TYPESENSE_HOST"),
#         "port": config("VITE_TYPESENSE_PORT"),
#         "protocol": "http"
#     }],
#     "api_key": config("VITE_TYPESENSE_KEY"),
#     "connection_timeout_seconds": 10
# }
#
# client = typesense.Client(TYPESENSE_CONFIG)
#
# # Retrieve schema of the collection
# schema = client.collections['autocomplete'].retrieve()
#
# # Get only the field names
# field_names = [field['name'] for field in schema['fields']]
#
# print(field_names)
#
# page = 1
# per_page = 100
# results = client.collections['autocomplete'].documents.search({
#     "q": "*",
#     "query_by": "id",  # Replace with any searchable field
#     "filter_by": "type:trait",
#     "per_page": per_page,
#     "page": page})
# hits = results.get("hits", [])
# print(hits)

#start = time.time()
#extract_variant_metrics(22, 16920281, "A", "T")
#end = time.time()
#print("It took ", end - start, " to retrieve all metrics for a single file!")