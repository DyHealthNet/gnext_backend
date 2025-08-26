import zarr
import numpy as np
import pandas as pd
import pysam
import re, os
from decouple import config
from django.conf import settings
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

from multiprocessing import Pool

logger = logging.getLogger("backend")


VARIANT_ARRAY = None

def init_worker(variant_array):
    global VARIANT_ARRAY
    VARIANT_ARRAY = variant_array

def get_all_variants_as_numpy():
    # Step 1: collect all variant IDs from your "master VCF"
    variants = []
    VCF_FILE = os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)
    with pysam.VariantFile(VCF_FILE) as vcf:
        for rec in vcf.fetch():
            vid = f"{rec.contig}:{rec.pos}:{rec.ref}:{rec.alts[0]}"
            variants.append(vid.encode())  # store as bytes for dtype="S"

    # Step 2: convert to NumPy fixed-width strings
    variant_array = np.array(variants, dtype="S40")  # 40 chars per variant is safe

    # Step 3: sort for searchsorted
    variant_array.sort()
    logger.info("Variant array constructed!")
    return variant_array

def process_trait(j, trait, filename):
    """
        Process one GWAS file and return the p-value column aligned to variant_array.
        """
    n = len(VARIANT_ARRAY)
    col_pval = np.full(n, np.nan, dtype="f4")

    with pysam.TabixFile(filename) as tbx:
        for k, line in enumerate(tbx.fetch()):
            f = line.strip().split("\t")
            # Example: CHR POS ID REF ALT ... PVAL (assume pval is in col 5)
            vid = f"{f[0]}:{f[1]}:{f[3]}:{f[4]}".encode()

            i = np.searchsorted(VARIANT_ARRAY, vid)
            if i < n and VARIANT_ARRAY[i] == vid:
                try:
                    pval = float(f[5])
                except Exception:
                    pval = np.nan
                col_pval[i] = pval

            if k > 0 and k % 1_000_000 == 0:
                logger.info(f"Worker trait={trait} processed {k:,} lines so far")

    return j, col_pval


def create_zarr_variants_trait_pvalues():
    # Get master variant index
    variant_array = get_all_variants_as_numpy()
    N_variants = len(variant_array)

    # Traits metadata
    pheno = pd.read_csv(config("PHENO_FILE"))
    traits = pheno["phenocode"].astype(str).tolist()
    paths = [
        os.path.join(
            settings.GWAS_NORM_DIR,
            re.sub(r'(\.[^.]+){1,2}$', '', os.path.basename(i)) + ".gz"
        )
        for i in pheno["filename"].astype(str).tolist()
    ]
    t2p = dict(zip(traits, paths))
    N_traits = len(traits)

    # Create Zarr store (row-chunked)
    root = zarr.open(settings.GWAS_ZARR_VARIANTS_TRAITS, mode="w")
    z_pvals = root.create_dataset(
        "pvalues",
        shape=(N_variants, N_traits),
        chunks=(10000, 1),  # row chunks (10k rows x 1 trait per chunk)
        dtype="f4",
        fill_value=np.nan,
        compressor=zarr.Blosc(cname="zstd", clevel=5, shuffle=2),
    )

    # Save metadata
    root.create_dataset("variants", data=variant_array)
    root.create_dataset("traits", data=np.array(traits, dtype="S40"))

    init_worker(variant_array)
    j, col_pval = process_trait(0, traits[0], t2p[traits[0]])
    logger.info("Non-NaN:", np.isfinite(col_pval).sum())

    # Run in parallel with ProcessPoolExecutor
    futures = {}
    with ProcessPoolExecutor(max_workers=int(config("MAX_WORKERS")),
                             initializer=init_worker,
                             initargs=(variant_array,)
                             ) as executor:
        for j, trait in enumerate(traits):
            fn = t2p[trait]
            future = executor.submit(process_trait, j, trait, fn)
            futures[future] = (j, trait)

        for future in as_completed(futures):
            j, trait = futures[future]
            col_pval = future.result()
            z_pvals[:, j] = col_pval
            logger.info(f"Stored trait {j + 1}/{N_traits}: {trait}")




