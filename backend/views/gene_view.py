from django.http import JsonResponse
from rest_framework import generics
import logging
import re

from backend.utils.extract_data_from_GWAS import extract_gene_signals, extract_region_associations
from backend.utils.typesense_client import get_gene_from_typesense

logger = logging.getLogger('backend')


class GeneOverview(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the PheWAS API.
        """
        gene_id = request.GET.get("id")
        gene_info = get_gene_from_typesense(gene_id)
        pattern = r"Chr\s*(\w+):\s*(\d+)\s*-\s*(\d+)\s*\(([-+])\)"
        m = re.match(pattern, gene_info["description"])
        if m:
            chrom, start, end, strand = m.groups()
            start = int(start)
            end = int(end)
        data = {
            "chromosome": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "symbol": gene_info["label"],
        }
        if data is None:
            return JsonResponse({"error": "No associations found for the given variant ID."}, status=404)
        else:
            return JsonResponse(data)

class GeneTopSignals(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests to the PheWAS API.
        """
        gene_id = request.GET.get("id")
        # Lookup gene positions
        gene_info = get_gene_from_typesense(gene_id)
        pattern = r"Chr\s*(\w+):\s*(\d+)\s*-\s*(\d+)\s*\(([-+])\)"
        m = re.match(pattern, gene_info["description"])
        if m:
            chrom, start, end, strand = m.groups()
            start = int(start)
            end = int(end)
        # Extract top signals
        top_signals = extract_gene_signals(chrom, start, end, strand, gene_id)
        if top_signals is None:
            return JsonResponse({"error": "No associations found for the given variant ID."}, status=404)
        else:
            data = {
                "rows": top_signals.to_dict(orient="records"),
                "header": top_signals.columns.tolist(),
            }
            return JsonResponse(data)

class GeneAssociationDataView(generics.GenericAPIView):
    def get(self, request, *args, **kwargs):
        """
        Handles GET requests for LocusZoom data.
        """
        trait_id = request.GET.get('trait_id')
        chr = request.GET.get('chr')
        start = request.GET.get('start')
        end = request.GET.get('end')
        
        if not all([trait_id, chr, start, end]):
            return JsonResponse({"error": "Missing required parameters"}, status=400)
        
        start = int(start)
        end = int(end)
        
        data = extract_region_associations(trait_id, chr, start, end)

        if data is None:
            return JsonResponse({"error": "No associations found for the given variant ID."}, status=404)
        else:
            return JsonResponse(data)