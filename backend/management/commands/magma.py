from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
import subprocess
from decouple import config
import gzip
from collections import defaultdict

logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting MAGMA execution of GWAS summary statistics files.")
           self.prepare_MAGMA_input()
           #self.run_MAGMA()
           logger.info("Finished MAGMA execution of GWAS summary statistics files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"MAGMA run with VEP failed: {e}")
           sys.exit(1)


    def extract_csq_fields(self, vcf_path):
        """
        Extracts the CSQ field names from a VCF file.
        """
        with gzip.open(vcf_path, 'rt') as f:
            for line in f:
                if line.startswith("##INFO") and "Format:" in line and "CSQ" in line:
                    parts = line.strip().split("|")
                    parts[0] = parts[0].split("Format: ")[1]
                    parts[-1] = parts[-1].strip('">')
                    return parts
        raise ValueError("No CSQ format found.")

    def prepare_MAGMA_input(self):
        GWAS_dir = config("GWAS_DIR")

        GWAS_vep_dir = os.path.join(GWAS_dir, "GWAS_vep")
        GWAS_annotated_vcf_file = os.path.join(GWAS_vep_dir, "annotated_full_variants.vcf.bgz")

        GWAS_magma_dir = os.path.join(GWAS_dir, "GWAS_magma")
        os.makedirs(GWAS_magma_dir, exist_ok=True)
        GWAS_magma_gene_annot_file = os.path.join(GWAS_magma_dir, "magma.genes.annot")

        without_gene_terms = ["regulatory_region_variant", "TF_binding_site_variant", "intergenic_variant",
                              "intron_variant"]

        csq_fields = self.extract_csq_fields(GWAS_annotated_vcf_file)

        window_size_up = config("MAGMA_WINDOW_UP")
        window_size_down = config("MAGMA_WINDOW_DOWN")

        gene_to_rsids = defaultdict(set)

        i = 1
        with gzip.open(GWAS_annotated_vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                info = parts[7]
                if "CSQ=" not in info:
                    continue
                csq_data = info.split("CSQ=")[1].split(";")[0]
                entries = csq_data.split(',')
                i += 1
                for entry in entries:
                    values = entry.split('|')
                    csq_dict = dict(zip(csq_fields, values))

                    if csq_dict.get("BIOTYPE") != "protein_coding" or csq_dict.get("Feature_type") != "Transcript":
                        continue

                    consequences = csq_dict.get("Consequence", "").split("&")
                    gene = csq_dict.get("Gene")
                    ids = csq_dict.get("Existing_variation", "").split("&")
                    rsid = next((i for i in ids if i.startswith("rs")), None)

                    if not gene or not rsid:
                        continue

                    # Check for consequences not in without_gene_terms
                    if any(c not in without_gene_terms for c in consequences):
                        gene_to_rsids[gene].add(rsid)
                        continue

                if i % 100000 == 0:
                    logger.info(f"Processed {i} lines from VCF file.")

        # Write MAGMA gene annotation file
        with open(GWAS_magma_gene_annot_file, 'w') as f:
            f.write("window_up = " + str(window_size_up) + "\n")
            for gene, rsids in gene_to_rsids.items():
                f.write(f"{gene}\t1:1:1\t{' '.join(rsids)}\n")
