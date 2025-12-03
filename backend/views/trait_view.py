from django.http import JsonResponse
from backend.utils.extract_data_from_GWAS import extract_variants_for_range, get_all_sign_variants_cutoff
from backend.utils.extract_data_from_manhattan import get_hits, wrap_generator_to_table_format
from backend.utils.converters import convert_variant_id
from rest_framework import generics
import logging
import re
import json
from django.conf import settings
from decouple import config
import os
import time
import pysam
import struct
import lmdb
import contextlib
import pandas as pd
import numpy as np

from backend.utils.typesense_client import get_phenotype_from_typesense

logger = logging.getLogger('backend')

class ManhattanView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the Manhattan API.
        """
        trait_id = request.GET.get("id")
        pheno_info = get_phenotype_from_typesense(trait_id)
        phenocode = pheno_info['id'] if pheno_info else None
        if phenocode is not None:
            manhattan_filepath = os.path.join(settings.MANHATTAN_DIR, phenocode + "_manhattan.json")
            with open(manhattan_filepath, "r") as f:
                manhattan_data = json.load(f)
            return JsonResponse(manhattan_data)
        else:
            logger.error(f"No phenotype found for trait: {trait_id}")
            return JsonResponse({"error": "Trait not found"}, status=404)

class QQView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the QQ API.
        """
        trait_id = request.GET.get("id")
        pheno_info = get_phenotype_from_typesense(trait_id)
        logger.info(pheno_info)
        phenocode = pheno_info['id'] if pheno_info else None
        if phenocode is not None:
            qq_filepath = os.path.join(settings.QQ_DIR, phenocode + "_qq.json")
            with open(qq_filepath, "r") as f:
                qq_data = json.load(f)
            return JsonResponse(qq_data)
        else:
            logger.error(f"No phenotype found for trait: {trait_id}")
            return JsonResponse({"error": "Trait not found"}, status=404)

class TraitInfoView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the QQ API.
        """
        trait_id = request.GET.get("id")
        pheno_info = get_phenotype_from_typesense(trait_id)
        logger.info(pheno_info)
        if pheno_info is not None:
            return JsonResponse(pheno_info)
        else:
            logger.error(f"No phenotype found for trait: {trait_id}")
            return JsonResponse({"error": "Trait not found"}, status=404)

class TraitView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to retrieve variant data for a given trait, chromosome, and position range.
        Supports queries by variant ID or chromosome range, with optional p-value cutoff filtering.
        """
        trait = request.GET.get("id")

        logger.info(f"Received trait request with trait: {trait}")
        pheno_info = get_phenotype_from_typesense(trait)
        phenocode = pheno_info['id'] if pheno_info else None
        if not phenocode:
            return JsonResponse({"error": "Trait not found"}, status=404)

        pval_cutoff_str = request.GET.get("pval_cutoff")
        if pval_cutoff_str in (None, "null", ""):
            pval_cutoff = 1.0
        else:
            pval_cutoff = float(request.GET.get("pval_cutoff", 1.0))

        mode = request.GET.get("mode")
        allowed_modes = {"loci", "pval", "rsid", "chromosome"}
        if mode not in allowed_modes:
            return JsonResponse({"error": "Invalid or missing 'mode' parameter. Must be one of: loci, pval, rsid, chromosome."}, status=400)
        if mode == "loci":
            generator = get_hits(pheno_info)
            data = wrap_generator_to_table_format(generator)
        elif mode == "pval":
            data = get_all_sign_variants_cutoff(phenocode, pval_cutoff=pval_cutoff)
        elif mode == "rsid":
            varid = request.GET.get("varid")
            if not varid:
                return JsonResponse({"error": "Missing required 'varid' parameter for rsid mode."}, status=400)
            neighbor_range = int(request.GET.get("range", 0))
            chr, pos, ref, alt = convert_variant_id(varid)
            start = max(pos - neighbor_range, 0)
            end = pos + neighbor_range
            data = extract_variants_for_range(phenocode, chr, start, end, pval_cutoff=pval_cutoff)
        elif mode == "chromosome":
            chr = request.GET.get("chr")
            if chr in (None, ""):
                return JsonResponse({"error": "Missing required 'chr' parameter for chromosome mode."}, status=400)
            start = int(request.GET.get("start", 0))
            end = int(request.GET.get("end", 0))
            data = extract_variants_for_range(phenocode, chr, start, end, pval_cutoff=pval_cutoff)

        if data is None:
            return JsonResponse({"error": "No variants found for the given trait."}, status=404)
        else:
            return JsonResponse(data)

class ChromosomeBoundsView(generics.GenericAPIView):

    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to retrieve the minimum and maximum position bounds for each chromosome
        in the GWAS file associated with a given trait. Returns a dictionary mapping chromosome names
        to their position bounds (min and max).
        """
        start_time = time.time()
        trait = request.GET.get("id")
        logger.info(f"Received chromosome request with trait: {trait}")
        pheno_info = get_phenotype_from_typesense(trait)
        phenocode = pheno_info['id'] if pheno_info else None
        if not phenocode:
            return JsonResponse({"error": "Trait not found"}, status=404)

        lmdb_path = settings.LMDB_RSID_FILE
        try:
            lmdb_env = lmdb.open(
                lmdb_path,
                map_size=1024 ** 4,
                max_dbs=25,
                readonly=True,
                lock=False,
                readahead=True,
                subdir=False
            )
            logger.debug("LMDB environment opened")
        except Exception as e:
            logger.warning(f"LMDB not available, chromosome bounds cannot be retrieved: {e}")
            lmdb_env = None

        bounds = {}

        # Open a read-only transaction
        with (lmdb_env.begin(buffers=True) if lmdb_env else contextlib.nullcontext()) as txn:
            # Iterate over chromosomes you care about
            for chrom in [str(i) for i in range(1, 23)] + ["X", "Y"]:
                try:
                    # Open DB handle for chromosome
                    db = lmdb_env.open_db(chrom.encode(), txn=txn)
                    cursor = txn.cursor(db)

                    # Get first key
                    if cursor.first():
                        min_pos = struct.unpack('>I', cursor.key())[0]

                        # Get last key
                        cursor.last()
                        max_pos = struct.unpack('>I', cursor.key())[0]

                        bounds[chrom] = {"min": min_pos, "max": max_pos}

                except lmdb.Error as e:
                    logger.warning(f"LMDB error for chromosome {chrom}: {e}")
                    # Skip chromosomes that don’t exist or other LMDB errors
                    continue

        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.debug(f"FINISHED FILTERING in {elapsed_time:.2f} seconds for getting the chromosome bounds")
        return JsonResponse(bounds, safe=False)

class MAGMAView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the MAGMA API.
        """
        trait_id = request.GET.get("id")
        if trait_id is not None:
            magma_filepath = os.path.join(settings.MAGMA_RESULTS_DIR, trait_id + "_magma.genes.out")
            magma_data = pd.read_csv(magma_filepath, delim_whitespace=True)
            # rename columns
            magma_data = magma_data.rename(columns={"GENE": "gene_id", "CHR": "chrom", "START": "start", "STOP": "end", "NSNPS": "#SNPs", "P": "pvalue", "ZSTAT": "zvalue", "SYMBOL":"gene_symbol"})
            # drop N and NPARAM column
            magma_data = magma_data.drop(columns=["N", "NPARAM"])
            # add -log10 pvalue column
            magma_data["-log10(pvalue)"] = -1 * np.log10(magma_data["pvalue"])
            magma_data["bonferroni_pvalue"] = magma_data["pvalue"] * len(magma_data)
            magma_data = magma_data.sort_values(by="pvalue")
            final_dict = {
                "rows": magma_data.to_dict(orient='records'),
                "header": list(magma_data.columns),
                "num_rows": len(magma_data)
            }
            return JsonResponse(final_dict)
        else:
            logger.error(f"No phenotype found for trait: {trait_id}")
            return JsonResponse({"error": "Trait not found"}, status=404)