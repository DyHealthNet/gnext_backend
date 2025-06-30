from django.http import JsonResponse
from decouple import config
from numpy import number
from rest_framework import generics
import logging
import re
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





