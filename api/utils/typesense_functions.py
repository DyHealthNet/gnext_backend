import typesense

def get_autocomplete_matches(query):
    client = typesense.Client({
        'api_key': 'xyz',
        'nodes': [{
            'host': 'localhost',
            'port': '8108',
            'protocol': 'http'
        }],
        'connection_timeout_seconds': 10
    })

    response = client.collections['phenotypes'].documents.search({
        'q': query,
        'query_by': 'label1',
        'per_page': 100
        })

    return(response)