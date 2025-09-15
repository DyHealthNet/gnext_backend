import time

from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import os
from decouple import config
import math
import json
import gzip
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed


from backend.utils.preprocessing.snp_to_rsid_mapping import setup_rsid_mapping_lmdb, add_rsID_with_lmdb
from backend.utils.preprocessing.locuszoom import manhattan, qq
from backend.utils.preprocessing.zorp.zorp import parsers, sniffers, readers, lookups
from backend.utils.preprocessing.magma.magma import read_magma_config
from django.conf import settings


logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting generation of Manhattan, QQ%s." % (" and MAGMA input files" if settings.MAGMA_ENABLED else " files"))
           self.generate_manhattan_qq_magma_files()
           if settings.MAGMA_ENABLED:
                  self.prepare_MAGMA_mapping_input()
           self.generate_top_hits()
           logger.info("Finished generation of Manhattan, QQ%s!" % (" and MAGMA input files" if settings.MAGMA_ENABLED else " files"))
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Generation of Manhattan{', QQ and MAGMA' if settings.MAGMA_ENABLED else ' and QQ'} files failed: {e}")
           sys.exit(1)

    @staticmethod
    def generate_manhattan_qq_magma_files():
        pheno_file = config("PHENO_FILE")

        GWAS_norm_dir = settings.GWAS_NORM_DIR
        GWAS_manhattan_dir = settings.GWAS_MANHATTAN_DIR
        os.makedirs(GWAS_manhattan_dir, exist_ok=True)
        GWAS_qq_dir = settings.GWAS_QQ_DIR
        os.makedirs(GWAS_qq_dir, exist_ok=True)
        GWAS_annotated_vcf_file = os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)
        lmdb_path = setup_rsid_mapping_lmdb(GWAS_annotated_vcf_file, GWAS_norm_dir)

        if settings.MAGMA_ENABLED:
            GWAS_magma_dir = settings.GWAS_MAGMA_DIR
            os.makedirs(GWAS_magma_dir, exist_ok=True)
            GWAS_magma_norm_dir = os.path.join(GWAS_magma_dir, "input_GWAS_norm")
            os.makedirs(GWAS_magma_norm_dir, exist_ok=True)
        else:
            GWAS_magma_norm_dir = None
            logger.debug("No need to prepare Magma input, since Magma is disabled.")

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        # Normalize GWAS files
        with ProcessPoolExecutor(max_workers=int(config("MAX_WORKERS"))) as executor:
            futures = [
                executor.submit(
                    Command.process_gwas_row,
                    r,
                    GWAS_manhattan_dir,
                    GWAS_qq_dir,
                    GWAS_norm_dir,
                    GWAS_magma_norm_dir,
                    lmdb_path
                )
                for i, r in pheno_dt.iterrows()
            ]
            for future in as_completed(futures):
                filename = future.result()
                logger.info(f"COMPLETED: Generated input file: %s", filename)

    @staticmethod
    def process_gwas_row(pheno, GWAS_manhattan_dir, GWAS_qq_dir, GWAS_norm_dir, GWAS_magma_norm_dir, lmdb_path):
        import logging
        logger = logging.getLogger("backend")

        norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(pheno['filename']))

        manhattan_filepath = os.path.join(GWAS_manhattan_dir, norm_filename + "_manhattan.json")
        qq_filepath = os.path.join(GWAS_qq_dir, norm_filename + "_qq.json")
        norm_filepath = os.path.join(GWAS_norm_dir, norm_filename + ".gz")
        # norm_with_rsid_filepath = norm_filepath.replace('.gz', '_rsid.gz')

        #TODO: Completed message also occurs when there was an error, change that

        # Update normalized GWAS files with rsID
        # if not os.path.exists(norm_with_rsid_filepath):
        #     logger.info("Started adding rsid to GWAS file: %s", norm_filepath)
        #     map_and_write_rsid(norm_filepath, lmdb_path)
        #     logger.info("COMPLETED: Added rsid to GWAS file: %s", norm_with_rsid_filepath)
        # else:
        #     logger.info("Skipping generation. GWAS file with rsid already exists: %s", norm_with_rsid_filepath)

        # Manhattan
        if not os.path.exists(manhattan_filepath):
            logger.info("Started Manhattan JSON file generation of GWAS file: %s", norm_filepath)
            reader_for_manhattan = sniffers.guess_gwas_standard(norm_filepath).add_filter('neg_log_pvalue')
            add_rsID_with_lmdb(reader_for_manhattan, lmdb_path)
            Command.generate_manhattan(reader_for_manhattan, manhattan_filepath)
            logger.info("COMPLETED: Manhattan JSON file generation of GWAS file: %s", norm_filepath)
        else:
            logger.info("Skipping generation. Manhattan JSON file already exists: %s", manhattan_filepath)

        # QQ
        if not os.path.exists(qq_filepath):
            logger.info("Startet QQ JSON file generation of GWAS file: %s", norm_filepath)
            reader_for_qq = sniffers.guess_gwas_standard(norm_filepath).add_filter('neg_log_pvalue')
            Command.generate_qq(reader_for_qq, qq_filepath)
            logger.info("COMPLETED: QQ JSON file generation of GWAS file: %s", norm_filepath)
        else:
            logger.info("Skipping generation. QQ JSON file already exists: %s", qq_filepath)

        # MAGMA
        if settings.MAGMA_ENABLED:
            magma_filepath = os.path.join(GWAS_magma_norm_dir, norm_filename + ".txt")
            if not os.path.exists(magma_filepath):
                logger.info("Started MAGMA normalized input file generation of GWAS file: %s", norm_filepath)
                reader_for_magma = sniffers.guess_gwas_standard(norm_filepath).add_filter('neg_log_pvalue')
                add_rsID_with_lmdb(reader_for_magma, lmdb_path)
                Command.generate_magma_input(reader_for_magma, magma_filepath, lmdb_path)
                logger.info("COMPLETED: MAGMA normalized input file generation of GWAS file: %s", norm_filepath)
            else:
                logger.info("Skipping generation. Magma input file already exists: %s", magma_filepath)

        return pheno['filename']

    @staticmethod
    def generate_manhattan(reader, out_filename: str) -> bool:
        """Generate manhattan plot data for the processed file"""

        binner = manhattan.Binner()
        for variant in reader:
            binner.process_variant(variant)

        manhattan_data = binner.get_result()

        # gl = get_genelocator(build, coding_only=False)
        for v_dict in manhattan_data['unbinned_variants']:
            # Annotate nearest gene(s) for all "top hits", and also clean up values so JS can handle them
            # It's possible to have more than one nearest gene for a given position (if variant is inside, not just near)
            # try:
            #    nearest_genes = [
            #        {
            #            'symbol': res['symbol'],
            #            'ensg': res['ensg']
            #       }
            #       for res in gl.at(v_dict["chrom"], v_dict["pos"])
            #   ]
            # except (gene_exc.BadCoordinateException, gene_exc.NoResultsFoundException):
            #   nearest_genes = []

            # v_dict['nearest_genes'] = nearest_genes

            if math.isinf(v_dict['neg_log_pvalue']):
                # JSON has no concept of infinity; use a string that browsers can type-coerce into the correct number
                v_dict['neg_log_pvalue'] = 'Infinity'

        with open(out_filename, 'w') as f:
            json.dump(manhattan_data, f)
        return True

    @staticmethod
    def generate_qq(reader, out_filename) -> bool:
        """Largely borrowed from PheWeb code (load.qq.make_json_file)"""
        # TODO: This step appears to load ALL data into memory (list on generator). This could be a memory hog; not sure if
        #   there is a way around it as it seems to rely on sorting values

        # TODO: Pheweb QQ code benefits from being passed { num_samples: n }, from metadata stored outside the
        #   gwas file. This is used when AF/MAF are present (which at the moment ingest pipeline does not support)
        variants = list(qq.augment_variants(reader))

        rv = {}
        if variants:
            if variants[0].maf is not None:
                rv['overall'] = qq.make_qq_unstratified(variants, include_qq=False)
                rv['by_maf'] = qq.make_qq_stratified(variants)
                rv['ci'] = list(qq.get_confidence_intervals(len(variants) / len(rv['by_maf'])))
            else:
                rv['overall'] = qq.make_qq_unstratified(variants, include_qq=True)
                rv['ci'] = list(qq.get_confidence_intervals(len(variants)))

        with open(out_filename, 'w') as f:
            json.dump(rv, f)
        return True

    @staticmethod
    def generate_magma_input(reader, out_filename, lmdb_path):
        start_time = time.time()
        reader.write(out_filename, columns=["rsid", "pval"], make_tabix=False)
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.debug(f"Time taken to normalize the GWAS stats file for Magma: {elapsed_time:.2f} seconds")

    @staticmethod
    def prepare_MAGMA_mapping_input():
        mconfig_rows = read_magma_config(config("MAGMA_CONFIG_FILE"))

        GWAS_annotated_vcf_file = settings.GWAS_ANNO_VCF_FILE
        GWAS_magma_dir = settings.GWAS_MAGMA_DIR
        GWAS_annotated_vcf_path = os.path.join(settings.GWAS_VEP_DIR, GWAS_annotated_vcf_file)

        for row in mconfig_rows:  # TODO use parallel processing
            mapping_strategy = (row.get("mapping_strategy") or "").lower()
            if mapping_strategy.strip() != "positional":
                return  # TODO adapt here depending on mapping strategy

            curr_window_up = row.get("window_up")
            curr_window_down = row.get("window_down")

            curr_GWAS_magma_dir = os.path.join(GWAS_magma_dir, f"magma_{mapping_strategy}_{curr_window_up}_{curr_window_down}")
            os.makedirs(curr_GWAS_magma_dir, exist_ok=True)

            GWAS_anno_magma_file = os.path.join(curr_GWAS_magma_dir, settings.GWAS_ANNO_MAGMA_FILE)

            # TODO if file of present has the correct window sizes in the header lines -> if not overwrite
            if not os.path.exists(GWAS_anno_magma_file):
                logger.info(f"Starting to create a SNP to gene annotation file for Magma...")

                without_gene_terms = ["regulatory_region_variant", "TF_binding_site_variant", "intergenic_variant",
                                      "intron_variant"]

                csq_fields = Command.extract_csq_fields(GWAS_annotated_vcf_path)

                gene_to_rsids = defaultdict(set)


                i = 1
                with gzip.open(GWAS_annotated_vcf_path, 'rt') as f:
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

                            if any(c not in without_gene_terms for c in consequences):
                                # If variant consequences not marked as up- or downstream it is inside the gene
                                if any(c not in ("upstream_gene_variant", "downstream_gene_variant") for c in
                                       consequences):
                                    gene_to_rsids[gene].add(rsid)
                                    continue
                                # If not check if variant falls in current window
                                distance = csq_dict.get("DISTANCE")
                                try:
                                    distance_int = int(distance)
                                except (ValueError, TypeError):
                                    continue
                                if "upstream_gene_variant" in consequences and curr_window_up != 0 and distance_int <= curr_window_up:
                                    gene_to_rsids[gene].add(rsid)
                                elif "downstream_gene_variant" in consequences and curr_window_down!=0 and distance_int <= curr_window_down:
                                    gene_to_rsids[gene].add(rsid)

                        if i % 100000 == 0:
                            logger.debug(f"Processed {i} lines from VCF file.")

                # Write MAGMA gene annotation file
                with open(GWAS_anno_magma_file, 'w') as f:
                    f.write("# window_up = " + str(curr_window_up) + "\n")
                    f.write("# window_down = " + str(curr_window_down) + "\n")
                    for gene, rsids in gene_to_rsids.items():
                        f.write(f"{gene}\t1:1:2\t{' '.join(rsids)}\n")

    @staticmethod
    def extract_csq_fields(vcf_path):
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

    @staticmethod
    def add_rsID_with_lmdb(reader, lmdb_path):
        build = os.path.join(lmdb_path, "data.mdb")
        rsid_finder = lookups.SnpToRsid(build, test=False)
        reader.add_lookup('rsid', lambda variant: rsid_finder(variant.chrom, variant.pos, variant.ref, variant.alt))

    @staticmethod
    def generate_top_hits(pval_cutoff = 1e-6, max_limit = 10000):
        """
        Generate a global top hits table across all phenotypes.
        Saves to GWAS_norm_dir/top_hits.json
        """
        pheno_file = config("PHENO_FILE")
        GWAS_norm_dir = settings.GWAS_NORM_DIR
        GWAS_manhattan_dir = settings.GWAS_MANHATTAN_DIR
        outpath = os.path.join(GWAS_norm_dir, "top_hits.json")

        if not os.path.exists(outpath):
            pheno_dt = pd.read_csv(pheno_file)
            phenos = pheno_dt.to_dict(orient="records")

            hits = []
            for pheno in phenos:
                hits.extend(get_hits(pheno, GWAS_manhattan_dir, pval_cutoff))

            # Sort by p-value ascending, keep top N
            hits.sort(key=lambda h: h["neg_log_pvalue"], reverse=True)
            hits = hits[:max_limit]

            # Keep only a subset of useful columns
            # cleaned_hits = [
            #     {k: hit.get(k) for k in ["category", "description", "pos", "pvalue", "phenocode", "top_variant", "category"]}
            #     for hit in hits
            # ]

            with open(outpath, "w") as f:
                json.dump(hits, f, indent=2)

            return hits

        return

def get_hits(pheno, GWAS_manhattan_dir, pval_cutoff):
    """
    Collect best variant per peak for a given phenotype.
    """
    norm_filename = re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(pheno["filename"]))
    manhattan_filepath = os.path.join(GWAS_manhattan_dir, norm_filename + "_manhattan.json")

    with open(manhattan_filepath) as f:
        variants = json.load(f)["unbinned_variants"]

    # Group by peak, keep lowest pvalue
    peak_to_best = {}
    for v in variants:
        if v.get("pvalue", 1.0) <= pval_cutoff and v.get("peak", False):
            # use chromosome+position as unique key
            key = (v["chrom"], v["pos"])
            best = peak_to_best.get(key)
            if best is None or v["pvalue"] < best["pvalue"]:
                v["phenocode"] = pheno["phenocode"]
                chrom = v["chrom"]
                pos = v["pos"]
                ref = v.get("ref", "")
                alt = v.get("alt", "")
                rsid = v.get("rsid")
                if rsid and rsid != ".":
                    v["top_variant"] = f"{chrom}_{pos}_{ref}/{alt} ({rsid})"
                else:
                    v["top_variant"] = f"{chrom}_{pos}_{ref}/{alt}"
                alt_allele_freq = v.get("alt_allele_freq")
                v["MAF"] = min(float(alt_allele_freq), 1 - float(alt_allele_freq)) if alt_allele_freq else None
                for k in ["description", "category"]:
                    if k in pheno:
                        v[k] = pheno[k]
                peak_to_best[key] = v

    yield from peak_to_best.values()
