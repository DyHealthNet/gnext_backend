from django.core.management.base import BaseCommand
import sys
import pandas as pd
import traceback
import logging
import os
from decouple import config
from backend.utils.preprocessing.locuszoom import manhattan, qq
import math
import json
from backend.utils.preprocessing.zorp.zorp import sniffers
from django.conf import settings


logger = logging.getLogger("backend")

class Command(BaseCommand):
    def handle(self, *args, **options):
       try:
           logger.info("Starting generation of Manhanttan and QQ files.")
           self.generate_manhattan_qq_files()
           logger.info("Finished generation of Manhattan and QQ files!")
       except Exception as e:
           # print stack trace
           traceback.print_exc()
           logger.error(f"Generation of Manhattan and QQ files failed: {e}")
           sys.exit(1)

    @staticmethod
    def generate_manhattan_qq_files():
        pheno_file = config("PHENO_FILE")

        GWAS_norm_dir = settings.GWAS_NORM_DIR
        GWAS_manhattan_dir = settings.GWAS_MANHATTAN_DIR
        os.makedirs(GWAS_manhattan_dir, exist_ok=True)
        GWAS_qq_dir = settings.GWAS_QQ_DIR
        os.makedirs(GWAS_qq_dir, exist_ok=True)

        # Importing phenotypes
        pheno_dt = pd.read_csv(pheno_file)

        # Normalize GWAS files
        for i, r in pheno_dt.iterrows():

            # Check if the Manhattan file already exists -> if yes, no need to process file again
            manhattan_filepath = os.path.join(GWAS_manhattan_dir, r['filename'].split(".")[0] + "_manhattan.json")
            qq_filepath = os.path.join(GWAS_qq_dir, r['filename'].split(".")[0] + "_qq.json")
            norm_filepath = os.path.join(GWAS_norm_dir, r['filename'] + ".gz")
            magma_filepath = "" # TODO: Bastienne

            if( os.path.exists(manhattan_filepath) and os.path.exists(qq_filepath)):
                logger.info("Skipping file %s, because input data (Manhattan, QQ, MAGMA) already present for file: ", r['filename'])
                continue
            else:
                # TODO: Bastienne -> add rsID to norm files

                # Generate Manhattan and QQ JSON file
                # Strong assumption: there are no invalid lines when a file reaches this stage; this operates on normalized data
                # Create two fresh readers
                reader_for_manhattan = sniffers.guess_gwas_standard(norm_filepath).add_filter('neg_log_pvalue')
                reader_for_qq = sniffers.guess_gwas_standard(norm_filepath).add_filter('neg_log_pvalue')

                Command.generate_manhattan(reader_for_manhattan, manhattan_filepath)
                logger.info("COMPLETED: Manhattan JSON file generation of GWAS file: %s", norm_filepath)
                Command.generate_qq(reader_for_qq, qq_filepath)
                logger.info("COMPLETED: QQ JSON file generation of GWAS file: %s", norm_filepath)

                # TODO: Bastienne -> add MAGMA input generation here

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

