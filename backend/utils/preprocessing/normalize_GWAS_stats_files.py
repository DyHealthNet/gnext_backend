# File for preprocessing GWAS summary statistics.
from backend.utils.preprocessing.zorp.zorp import sniffers
from backend.utils.preprocessing.zorp.zorp import parsers


def normalize_gwas_stats_file(input_file: str, output_file: str) -> None:
    """
    Preprocess GWAS summary statistics file to normalize and filter data.

    Args:
        input_file (str): Path to the input GWAS summary statistics file.
        output_file (str): Path to save the preprocessed output file.
    """

    parser_options = { #TODO: make these configurable
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

    reader = sniffers.guess.gwas_generic(input_file, parser = parser, skip_errors = True)
    reader.write(output_file, make_tabix=True)

