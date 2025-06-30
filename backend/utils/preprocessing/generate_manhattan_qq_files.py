#Extracted from locuszoom-hosted: d726fb2
from backend.utils.preprocessing.zorp.zorp import sniffers
import json
import math
from backend.utils import preprocessing as manhattan, preprocessing as qq


def generate_manhattan(reader, out_filename: str) -> bool:
    """Generate manhattan plot data for the processed file"""

    binner = manhattan.Binner()
    for variant in reader:
        binner.process_variant(variant)

    manhattan_data = binner.get_result()

    #gl = get_genelocator(build, coding_only=False)
    for v_dict in manhattan_data['unbinned_variants']:
        # Annotate nearest gene(s) for all "top hits", and also clean up values so JS can handle them
        # It's possible to have more than one nearest gene for a given position (if variant is inside, not just near)
        #try:
        #    nearest_genes = [
        #        {
        #            'symbol': res['symbol'],
        #            'ensg': res['ensg']
         #       }
         #       for res in gl.at(v_dict["chrom"], v_dict["pos"])
         #   ]
        #except (gene_exc.BadCoordinateException, gene_exc.NoResultsFoundException):
         #   nearest_genes = []

        #v_dict['nearest_genes'] = nearest_genes

        if math.isinf(v_dict['neg_log_pvalue']):
            # JSON has no concept of infinity; use a string that browsers can type-coerce into the correct number
            v_dict['neg_log_pvalue'] = 'Infinity'

    with open(out_filename, 'w') as f:
        json.dump(manhattan_data, f)
    return True

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

def generate_manhattan_qq_json(in_filename: str, manhattan_file: str, qq_file: str) -> bool:
    # Strong assumption: there are no invalid lines when a file reaches this stage; this operates on normalized data
    reader = sniffers.guess_gwas_standard(in_filename) \
        .add_filter('neg_log_pvalue')

    generate_manhattan(reader, manhattan_file)
    generate_qq(reader, qq_file)
    return True
