import typesense
import gzip
import json
import csv
import pandas as pd

from decouple import config

# Collection name
collection_name = "entries"

# Config variables
api_key = config('TYPESENSE_KEY')
typesense_host = config('TYPESENSE_HOST')
typesense_port = config('TYPESENSE_PORT')

chr_column = config('CHR_COLUMN')
pos_column = config('POS_COLUMN')
ref_column = config('REF_COLUMN')
alt_column = config('ALT_COLUMN')

client = typesense.Client({
    'api_key': api_key,
    'nodes': [{
        'host': typesense_host,
        'port': typesense_port,
        'protocol': 'http'
    }],
    'connection_timeout_seconds': 10
})


# count number of entries in variants collection
print(client.collections.retrieve())
