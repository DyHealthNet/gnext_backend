from typing import Dict, List, Any, Mapping, Sequence

from django.conf import settings
import typesense
import logging


logger = logging.getLogger('backend')

client = typesense.Client(settings.TYPESENSE_CONFIG)

def get_client():
    """
    Returns the Typesense client.
    """
    return client

def get_phenotype_from_typesense(trait_id):
    client = get_client()
    doc = client.collections['autocomplete'].documents[trait_id].retrieve()
    return(doc)


def get_all_phenotypes_from_typesense():
    client = get_client()

    all_documents = []
    page = 1
    per_page = 100

    try:
        while True:
            results = client.collections['autocomplete'].documents.search({
                "q": "*",
                "query_by": "description",  # Replace with any searchable field
                "filter_by": "type:trait",
                "per_page": per_page,
                "page": page
            })

            hits = results.get("hits", [])
            if not hits:
                break

            all_documents.extend([hit["document"] for hit in hits])
            page += 1

        return all_documents

    except Exception as e:
        return None

def get_number_of_documents():
    client = get_client()
    results = client.collections['autocomplete'].documents.search({
        "q": "*",
        "query_by": "description",
        "facet_by": "type",
        "per_page": 0,

    })
    type_counts = {entry['value']: entry['count'] for entry in results['facet_counts'][0]['counts']}
    return type_counts