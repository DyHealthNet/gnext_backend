from django.http import JsonResponse, HttpResponse
from decouple import config
from numpy import number
from rest_framework import generics
import logging
import re
import json
from django.conf import settings
import os
from backend.utils.typesense_client import get_number_of_documents

logger = logging.getLogger('backend')

class OverviewDataView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        number_docs = get_number_of_documents()
        if number_docs is None:
            return JsonResponse({"error": "Failed to retrieve data"}, status=500)
        else:
            # Return the number of documents as a JSON response
            return JsonResponse(number_docs, status=200)


class TopHitsView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to retrieve top hits data from a JSON file and returns it as a JSON response.
        """
        logger.debug(f"top_hits_file: {settings.TOP_HITS_FILE}")
        try:
            with open(settings.TOP_HITS_FILE, "r") as f:
                top_hits = json.load(f)
                # Make neg_log_pvalue = Infinity to string for JSON serialization
                json_str = json.dumps(top_hits).replace('Infinity', '"Infinity"')
            return HttpResponse(json_str, content_type="application/json")
        except Exception as e:
            logger.error(f"Error opening Top Hits file {settings.TOP_HITS_FILE}: {e}")
            return JsonResponse({"error": "Failed to open Top Hits file"}, status=500)

class ConfigView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the MAGMA API.
        """
        try:
            # Read JSON file of Nextflow parameters -> NF_PARAM_FILE
            with open(settings.NF_PARAM_FILE, "r") as f:
                nf_params = json.load(f)
            magma_ref_pop = config("MAGMA_REF_POPULATION")
            magma_ref_gene_loc = config("MAGMA_REF_GENE_LOCATION")
            nf_params["magma_ref_pop"] = magma_ref_pop
            nf_params["magma_ref_gene_loc"] = magma_ref_gene_loc
            return JsonResponse(nf_params)
        except:
            return JsonResponse({"error": "Failed to retrieve data"}, status=500)

