from django.http import JsonResponse

from backend.utils.converters import convert_variant_id
from backend.utils.extract_data_from_GWAS import extract_variant_metrics
from backend.utils.extract_data_from_VEP import extract_variant_annotation
from decouple import config
from rest_framework import generics
import logging
import re

logger = logging.getLogger('backend')

class VariantMetricsView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the PheWAS API.
        """
        build = request.GET.get("build")
        filter = request.GET.get("filter")
        format = request.GET.get("format")

        logger.info(f"Received request with build: {build}, filter: {filter}, format: {format}")

        variant_id = re.search(r"'([^']+)'", filter).group(1)
        chr, pos, ref, alt = convert_variant_id(variant_id)
        results, min_af, max_af = extract_variant_metrics(chr, pos, ref, alt)
        json_resp = {"data": results,
                     "lastPage": None,
                     "meta": {
                         "build": [config("VITE_GENOME_BUILD")],
                         "min_af": min_af if min_af != float("inf") else None,
                         "max_af": max_af if max_af != float("-inf") else None}}
        if results is None:
            return JsonResponse({"error": "No associations found for the given variant ID."}, status=404)
        else:
            return JsonResponse(json_resp)

class VariantAnnotationView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):

        variant_id = request.GET.get("id")

        logger.info(f"Received request with id for variant annotation: {variant_id}")

        data = extract_variant_annotation(variant_id)
        if data is None:
            return JsonResponse({"error": "No annotation found for the given variant ID."}, status=404)
        else:
            return JsonResponse(data)