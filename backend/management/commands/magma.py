from django.core.management.base import BaseCommand
import sys
import traceback
import logging
import os
import time
import subprocess
from decouple import config
import gzip
from collections import defaultdict
import pandas as pd
from backend.utils.preprocessing.zorp.zorp import parsers, sniffers
import backend.utils.preprocessing.magma.magma_norm_exec as magma_exec

logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting MAGMA execution of GWAS summary statistics files.")
           self.prepare_MAGMA_input()
           self.normalze_MAGMA_input()
           self.run_MAGMA()
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
        GWAS_annotated_vcf_file = os.path.join(GWAS_vep_dir, "annotated_variants.vcf.bgz") # annotated_full_variants

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
            f.write("# window_up = " + str(window_size_up) + "\n")
            f.write("# window_down = " + str(window_size_down) + "\n")
            for gene, rsids in gene_to_rsids.items():
                f.write(f"{gene}\t1:1:2\t{' '.join(rsids)}\n")

    def normalze_MAGMA_input(self):
        dir_path = config('GWAS_DIR')
        vep_path = config('VEP_ANNO_FILE')
        pheno_file = config('PHENO_FILE')
        pheno_dt = pd.read_csv(pheno_file)

        parser_options = {
            "chrom_col": 1,
            "pos_col": 2,
            "ref_col": 3,
            "alt_col": 4,
            "pval_col": 8,
            "is_neg_log_pvalue": True,
            'beta': 6,
            'stderr_beta': 7,
            'alt_allele_freq': 5,
            'rsid': None
        }

        parser = parsers.GenericGwasLineParser(**parser_options)

        os.makedirs(dir_path + "GWAS_stats_norm/", exist_ok=True)
        os.makedirs(dir_path + "GWAS_magma/", exist_ok=True)

        lmdb_path = dir_path + "/GWAS_magma/" + "lmdb_" + config('VITE_GENOME_BUILD')
        # Only build the LMDB if it doesn't exist or is missing required files
        if not os.path.isdir(lmdb_path) or not os.path.exists(os.path.join(lmdb_path, "data.mdb")):
            logger.info("LMDB not found, creating...")
            start_time = time.time()
            magma_exec.build_snp_map_lmdb_from_vcf(vep_path, lmdb_path, map_size=10 * 1024 ** 3)
            end_time = time.time()
            logger.debug(f"Time taken to produce LMDB mapping lib: {end_time - start_time:.2f} seconds")
        else:
            logger.info("LMDB already exists, skipping creation.")
        # Time taken to produce Lmdb mapping lib: 529.44 seconds

        # Importing phenotypes
        for i, r in pheno_dt.iterrows():
            logger.debug(f"Importing phenotype: {r['phenocode']}")
            GWAS_file = r['filename']
            sample_file = dir_path + GWAS_file
            norm_gwas_file = dir_path + "/GWAS_magma/GWAS_stats_norm/" + GWAS_file.replace(".tsv.bgz", ".txt")
            if not os.path.exists(norm_gwas_file):

                reader = sniffers.guess_gwas_generic(sample_file, parser=parser, skip_errors=True)

                start_time = time.time()
                status = magma_exec.normalize_contents_lib(reader,
                                                           norm_gwas_file, genome_build='GRCh37',
                                                           debug_mode=False, lmdb_path=lmdb_path)
                end_time = time.time()
                elapsed_time = end_time - start_time
                logger.debug(f"Time taken to normalize the GWAS stats file: {elapsed_time:.2f} seconds")
            else:
                logger.debug(f"Normalized GWAS stats file already present: {GWAS_file}")

    def run_MAGMA(self):
        magma = config('MAGMA_EXEC')
        dir_path = config('GWAS_DIR')
        gene_annot = os.path.join(dir_path, "GWAS_magma/magma.genes.annot")
        LD_path = config('MAGMA_LD_REF')
        pheno_file = config('PHENO_FILE')
        pheno_dt = pd.read_csv(pheno_file)

        n_samples = config('N_SAMPLES')
        for i, r in pheno_dt.iterrows():
            gwas_file = r['filename']
            sample_file = dir_path  + "/GWAS_magma/GWAS_stats_norm/" + gwas_file.replace(".tsv.bgz", ".txt")
            magma_file = dir_path + "/GWAS_magma/MAGMA_results/" + gwas_file.replace(".tsv.bgz", "")

            if not os.path.exists(magma_file):

                gene_cmd = [
                    magma,
                    '--bfile', LD_path,
                    'synonym-dup=drop-dup',  # optional or change!!
                    '--gene-annot', gene_annot,
                    '--pval', sample_file, f'N={n_samples}',
                    '--gene-model', 'snp-wise=mean',  # optional or change but this ones best and fastest for summary stats
                    '--genes-only',  # since we do not perform pathway analysis with magma we don't need this
                    '--out', magma_file
                ]

                logger.info(f"[INFO] Running gene analysis for {gwas_file}")
                start_time = time.time()
                result = subprocess.run(gene_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                logger.debug("MAGMA stdout:\n%s", result.stdout)
                logger.debug("MAGMA stderr:\n%s", result.stderr)
                if result.returncode != 0:
                    logger.error("MAGMA failed with exit code %d", result.returncode)
                    raise subprocess.CalledProcessError(result.returncode, gene_cmd, output=result.stdout,
                                                        stderr=result.stderr)
                end_time = time.time()
                elapsed_time = end_time - start_time
                logger.debug(f"[DONE] Finished MAGMA run {magma_file} in {elapsed_time:.2f} seconds")