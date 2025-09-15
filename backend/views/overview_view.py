from django.http import JsonResponse
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
        GWAS_norm_dir = settings.GWAS_NORM_DIR
        top_hits_file = os.path.join(GWAS_norm_dir, "top_hits.json")
        logger.debug(f"top_hits_file: {top_hits_file}")
        try:
            with open(top_hits_file, "r") as f:
                top_hits = json.load(f)
            return JsonResponse(top_hits, safe=False)
        except Exception as e:
            logger.error(f"Error opening Top Hits file {top_hits_file}: {e}")
            return JsonResponse({"error": "Failed to open Top Hits file"}, status=500)





