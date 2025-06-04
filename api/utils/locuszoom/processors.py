#Extracted from locuszoom-hosted: d726fb2
from api.utils.zorp.zorp import parsers, sniffers, readers, lookups
import json
import math
import api.utils.locuszoom.manhattan as manhattan
import api.utils.locuszoom.qq as qq
from genelocator import get_genelocator
import genelocator.exception as gene_exc


def normalize_contents(reader: readers.BaseReader, dest_path: str, build: str, debug_mode=False) -> bool:
    """
    Initial content ingestion: load the file and write variants in a standardized format

    This routine will deliberately exclude lines that could not be handled in a reliable fashion, such as pval=NA

    In "debug mode", the ingest process will use a smaller (test environment optimized) version of the
        rsid_finder lookup.
    """
    #rsid_finder = lookups.SnpToRsid(build, test=debug_mode)
    #reader.add_lookup('rsid', lambda variant: rsid_finder(variant.chrom, variant.pos, variant.ref, variant.alt))
    reader.write(dest_path, make_tabix=True)
    # In reality a failing task will usually raise an exception rather than returning False
    return True

def generate_manhattan(build: str, in_filename: str, out_filename: str) -> bool:
    """Generate manhattan plot data for the processed file"""
    # Strong assumption: there are no invalid lines when a file reaches this stage; this operates on normalized data
    reader = sniffers.guess_gwas_standard(in_filename)\
        .add_filter('neg_log_pvalue')

    binner = manhattan.Binner()
    for variant in reader:
        binner.process_variant(variant)

    manhattan_data = binner.get_result()

    gl = get_genelocator(build, coding_only=False)
    for v_dict in manhattan_data['unbinned_variants']:
        # Annotate nearest gene(s) for all "top hits", and also clean up values so JS can handle them
        # It's possible to have more than one nearest gene for a given position (if variant is inside, not just near)
        try:
            nearest_genes = [
                {
                    'symbol': res['symbol'],
                    'ensg': res['ensg']
                }
                for res in gl.at(v_dict["chrom"], v_dict["pos"])
            ]
        except (gene_exc.BadCoordinateException, gene_exc.NoResultsFoundException):
            nearest_genes = []

        v_dict['nearest_genes'] = nearest_genes

        if math.isinf(v_dict['neg_log_pvalue']):
            # JSON has no concept of infinity; use a string that browsers can type-coerce into the correct number
            v_dict['neg_log_pvalue'] = 'Infinity'

    with open(out_filename, 'w') as f:
        json.dump(manhattan_data, f)
    return True

def generate_qq(in_filename: str, out_filename) -> bool:
    """Largely borrowed from PheWeb code (load.qq.make_json_file)"""
    # TODO: This step appears to load ALL data into memory (list on generator). This could be a memory hog; not sure if
    #   there is a way around it as it seems to rely on sorting values
    reader = sniffers.guess_gwas_standard(in_filename)\
        .add_filter("neg_log_pvalue")

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

