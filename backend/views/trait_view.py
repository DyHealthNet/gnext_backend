from django.http import JsonResponse
from backend.utils.extract_data_from_GWAS import extract_phenotype_results_for_variant
from backend.utils.extract_data_from_VEP import extract_variant_annotation
from decouple import config
from rest_framework import generics
import logging
import re
import json
from django.conf import settings
import os


from backend.utils.typesense_client import get_phenotype_from_typesense

logger = logging.getLogger('backend')

class ManhattanView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the Manhattan API.
        """
        trait = request.GET.get("trait")
        dir_path = config("GWAS_DIR")
        logger.info(f"Received request with trait: {trait}")
        pheno_info = get_phenotype_from_typesense(trait)
        file_name = pheno_info[0]['filename'] if pheno_info else None
        if file_name is not None:
            norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(file_name))
            manhattan_filepath = os.path.join(settings.GWAS_MANHATTAN_DIR, norm_filename + "_manhattan.json")
            with open(manhattan_filepath, "r") as f:
                manhattan_data = json.load(f)
            return JsonResponse(manhattan_data)
        else:
            logger.error(f"No phenotype found for trait: {trait}")
            return JsonResponse({"error": "Trait not found"}, status=404)

class QQView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the QQ API.
        """
        trait = request.GET.get("trait")
        dir_path = config("GWAS_DIR")
        logger.info(f"Received request with trait: {trait}")
        pheno_info = get_phenotype_from_typesense(trait)
        file_name = pheno_info[0]['filename'] if pheno_info else None
        if file_name is not None:
            norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(file_name))
            qq_filepath = os.path.join(settings.GWAS_QQ_DIR, norm_filename + "_qq.json")
            with open(qq_filepath, "r") as f:
                qq_data = json.load(f)
            return JsonResponse(qq_data)
        else:
            logger.error(f"No phenotype found for trait: {trait}")
            return JsonResponse({"error": "Trait not found"}, status=404)