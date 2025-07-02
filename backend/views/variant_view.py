from django.http import JsonResponse
from backend.utils.extract_data_from_GWAS import extract_phenotype_results_for_variant
from backend.utils.extract_data_from_VEP import extract_variant_annotation
from decouple import config
from rest_framework import generics
import logging
import re

logger = logging.getLogger('backend')

class PheWASView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the PheWAS API.
        """
        build = request.GET.get("build")
        filter = request.GET.get("filter")
        format = request.GET.get("format")

        variant_id = re.search(r"'([^']+)'", filter).group(1)

        logger.info(f"Received request with build: {build}, filter: {filter}, format: {format}")

        data = extract_phenotype_results_for_variant(variant_id)
        if data is None:
            return JsonResponse({"error": "No associations found for the given variant ID."}, status=404)
        else:
            return JsonResponse({
                "data": data,
                "lastPage": None,
                "meta": {
                    "build": [
                        "GRCh37"  # TODO: extract from configuration file
                    ]
                }
            })

class VariantAnnotationView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):

        variant_id = request.GET.get("id")

        logger.info(f"Received request with id for variant annotation: {variant_id}")

        data = extract_variant_annotation(variant_id)
        if data is None:
            return JsonResponse({"error": "No annotation found for the given variant ID."}, status=404)
        else:
            return JsonResponse(data)