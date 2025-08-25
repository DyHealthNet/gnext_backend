import pysam
from django.http import JsonResponse
from backend.utils.extract_data_from_GWAS import extract_phenotype_results_for_variant, extract_variants_for_range
from backend.utils.extract_data_from_VEP import extract_variant_annotation
from backend.utils.converters import convert_variant_id
from decouple import config
from rest_framework import generics
import logging
import re
import json
from django.conf import settings
import os

import time


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

class TraitView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to retrieve variant data for a given trait, chromosome, and position range.
        Supports queries by variant ID or chromosome range, with optional p-value cutoff filtering.
        """
        trait = request.GET.get("trait")
        # Get query parameters
        pval_cutoff_str = request.GET.get("pval_cutoff")
        if pval_cutoff_str in (None, "null", ""):
            pval_cutoff = 1.0
        else:
            pval_cutoff = float(request.GET.get("pval_cutoff", 1.0))

        # rsid mode
        varid = request.GET.get("varid")
        neighbor_range = int(request.GET.get("range", 0))

        if varid:
            chr, pos, ref, alt = convert_variant_id(varid)
            start = max(pos - neighbor_range, 0)
            end = pos + neighbor_range
        else:
            # chromosome range mode
            chr = request.GET.get("chr")
            start = int(request.GET.get("start",0))
            end = int(request.GET.get("end",0))

        logger.info(f"Received request with trait: {trait}")
        pheno_info = get_phenotype_from_typesense(trait)
        file_name = pheno_info[0]['filename'] if pheno_info else None
        if not file_name:
            return JsonResponse({"error": "Trait not found"}, status=404)

        data = extract_variants_for_range(file_name, chr, start, end, pval_cutoff=pval_cutoff)

        if data is None:
            return JsonResponse({"error": "No variants found for the given trait."}, status=404)
        else:
            return JsonResponse(data)

class ChromosomeBoundsView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to retrieve variant data for a given trait, chromosome, and position range.
        Supports queries by variant ID or chromosome range, with optional p-value cutoff filtering.
        """
        start_time = time.time()
        trait = request.GET.get("trait")
        logger.info(f"Received request with trait: {trait}")
        pheno_info = get_phenotype_from_typesense(trait)
        file_name = pheno_info[0]['filename'] if pheno_info else None
        if not file_name:
            return JsonResponse({"error": "Trait not found"}, status=404)

        norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(file_name))
        norm_filepath = os.path.join(settings.GWAS_NORM_DIR, norm_filename + ".gz")

        try:
            tabix_file = pysam.TabixFile(norm_filepath)
        except Exception as e:
            logger.error(f"Error opening Tabix file {norm_filepath}: {e}")
            return JsonResponse({"error": "Failed to open GWAS file"}, status=500)

        bounds = {}
        for chrom in tabix_file.contigs:
            try:
                # Fetch the first record
                first = next(tabix_file.fetch(chrom))
                min_pos = int(first.split("\t")[1])

                # Fetch the last record by iterating from the end (pysam doesn’t have direct last fetch, but can use reversed iterator)
                last = None
                for last in tabix_file.fetch(chrom):
                    pass
                max_pos = int(last.split("\t")[1]) if last else min_pos

                bounds[chrom] = {"min": min_pos, "max": max_pos}
            except (StopIteration, ValueError):
                # chromosome has no records
                continue

        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.debug(f"FINISHED FILTERING in {elapsed_time:.2f} seconds for getting the chromosome bounds")

        return JsonResponse(bounds, safe=False)