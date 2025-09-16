import os
import json
import logging
import re

logger = logging.getLogger("backend")

def get_hits(pheno, GWAS_manhattan_dir, pval_cutoff = 1e-6):
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
                # typesense turns key phenocode into id?
                if "phenocode" in pheno:
                    v["phenocode"] = pheno["phenocode"]
                elif "id" in pheno:
                    v["id"] = pheno["id"]
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

def wrap_generator_to_table_format(generator):
    variants = list(generator)

    if not variants:
        return {"rows": {}, "header": []}

    # Define the columns you want in the table (and in order)
    desired_headers = [
        "top_variant",
        "neg_log_pvalue",
        "beta",
        "stderr_beta",
        "alt_allele_freq",
        "pvalue",
        "peak",
    ]

    # filter each variant down to only the desired keys
    filtered_variants = [
        {k: v.get(k) for k in desired_headers if k in v}
        for v in variants
    ]

    rows = {i: v for i, v in enumerate(filtered_variants)}
    return {"rows": rows, "header": desired_headers}