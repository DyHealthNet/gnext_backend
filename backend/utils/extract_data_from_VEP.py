from decouple import config
import pysam
from backend.utils.converters import convert_variant_id
import gzip
import logging
import pandas as pd
import os
from backend.utils.VEP_consequences import VEP_RANK_DICT
from django.conf import settings

logger = logging.getLogger('backend')

def get_most_severe(consequences):
    """
        Get the most severe consequence term from a list.
        Returns None if no recognized terms are found.
        """
    most_severe = None
    most_severe_impact = None
    best_rank = float('inf')

    for term in consequences:
        if "&" in term:
            terms = term.split("&")
        else:
            terms = [term]
        for t in terms:
            rank = VEP_RANK_DICT[t]["rank"]
            if rank is not None and rank < best_rank:
                best_rank = rank
                most_severe = t
                most_severe_impact = VEP_RANK_DICT[t]["impact"]
    return most_severe, most_severe_impact


def extract_variant_annotation(variant_id):
    vep_anno_file = os.path.join(settings.GWAS_VEP_DIR, settings.GWAS_ANNO_VCF_FILE)

    # Extract GWAS results of variant from phenotype file via tabix
    tabix_file = pysam.TabixFile(vep_anno_file)
    # Get variant information
    logger.info("INFO: Extracting variant annotation for variant ID: " + variant_id)
    chr, pos, ref, alt = convert_variant_id(variant_id)

    with gzip.open(vep_anno_file, 'rt') as vcf:
        for line in vcf:
            if line.startswith("##INFO"):
                anno_columns = line.strip().split("|")
                anno_columns[0] = anno_columns[0].split("Format: ")[1].split('"')[0]
                anno_columns[len(anno_columns) - 1] = anno_columns[len(anno_columns) - 1].strip('">')
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                columns = [h.replace("#", "") for h in header]
                break

    try:
        allele_frequencies = {}
        external_ids = []
        rows = []
        consequences = []
        location = chr + ":" + str(pos)
        AF = 0
        closest_gene = []
        transcript_dict = {}
        regulatory_dict = {}
        motif_dict = {}
        for row in tabix_file.fetch(chr, pos - 1, pos):
            row = row.split("\t")
            row = dict(zip(columns, row))
            info = row["INFO"].replace("CSQ=", "").split(",")
            info_dict = [dict(zip(anno_columns, i.split("|"))) for i in info]
            info_pd = pd.DataFrame(info_dict)

            # Allele frequencies
            af_keys = [key for key in anno_columns if key.endswith("_AF") and key != "MAX_AF"]
            allele_frequencies = info_pd[["AF"] + af_keys].replace('', pd.NA).dropna(axis=1, how='all').astype('float').to_dict(orient='records')[0]
            # rename keys to be more readable
            new_af_names = {"AF": "Global AF (1000 Genomes)",
                            "AFR_AF": "AFR (1000 Genomes)",
                            "AMR_AF": "AMR (1000 Genomes)",
                            "EAS_AF": "EAS (1000 Genomes)",
                            "EUR_AF": "EUR (1000 Genomes)",
                            "SAS_AF": "SAS (1000 Genomes)",
                            "gnomADe_AF": "Global AF (gnomAD Exome)",
                            "gnomADe_AFR_AF": "AFR (gnomAD Exome)",
                            "gnomADe_AMR_AF": "AMR (gnomAD Exome)",
                            "gnomADe_ASJ_AF": "ASJ (gnomAD Exome)",
                            "gnomADe_EAS_AF": "EAS (gnomAD Exome)",
                            "gnomADe_FIN_AF": "FIN (gnomAD Exome)",
                            "gnomADe_NFE_AF": "NFE (gnomAD Exome)",
                            "gnomADe_OTH_AF": "OTH (gnomAD Exome)",
                            "gnomADe_SAS_AF": "SAS (gnomAD Exome)",
                            "gnomADg_AF": "Global AF (gnomAD Genome)",
                            "gnomADg_AFR_AF": "AFR (gnomAD Genome)",
                            "gnomADg_AMR_AF": "AMR (gnomAD Genome)",
                            "gnomADg_ASJ_AF": "ASJ (gnomAD Genome)",
                            "gnomADg_EAS_AF": "EAS (gnomAD Genome)",
                            "gnomADg_FIN_AF": "FIN (gnomAD Genome)",
                            "gnomADg_NFE_AF": "NFE (gnomAD Genome)",
                            "gnomADg_OTH_AF": "OTH (gnomAD Genome)",
                            "gnomADg_SAS_AF": "SAS (gnomAD Genome)",
                            "gnomADg_MID_AF": "MID (gnomAD Genome)",
                            "gnomADg_AMI_AF": "AMI (gnomAD Genome)",
                            "gnomADg_REMAINING_AF": "Remaining (gnomAD Genome)",
                            }

            allele_frequencies = {new_af_names.get(k, k): v for k, v in allele_frequencies.items()}


            external_ids = info_pd["Existing_variation"].str.split("&").explode().drop_duplicates().tolist()
            consequences = info_pd["Consequence"].drop_duplicates().tolist()
            # transcript consequences
            transcript_columns = ["Feature", "BIOTYPE", "Gene", "SYMBOL", "IMPACT", "Consequence", "DISTANCE","ENSP", "EXON", "INTRON", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids"]

            transcript_consequences_dt = info_pd[info_pd["Feature_type"] == "Transcript"][transcript_columns]
            transcript_consequences_dt["Consequence"] = transcript_consequences_dt["Consequence"].apply(
                lambda x: x.replace("&", ", ") if isinstance(x, str) else x
            )
            transcript_dict = {"headers": transcript_columns, "rows": transcript_consequences_dt.to_dict(orient="records")}

            # regulatory consequences
            regulatory_columns = ["Feature", "BIOTYPE", "Consequence"]
            regulatory_consequences_dt = info_pd[info_pd["Feature_type"] == "RegulatoryFeature"][regulatory_columns]
            regulatory_consequences_dt["Consequence"] = regulatory_consequences_dt["Consequence"].apply(
                lambda x: x.replace("&", ", ") if isinstance(x, str) else x
            )
            regulatory_dict = {"headers": regulatory_columns, "rows": regulatory_consequences_dt.to_dict(orient="records")}

            # motif consequences
            motif_columns = ["Feature", "Consequence", "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS"]
            motif_consequences_dt = info_pd[info_pd["Feature_type"] == "MotifFeature"][motif_columns]
            motif_consequences_dt["Consequence"] = motif_consequences_dt["Consequence"].apply(
                lambda x: x.replace("&", ", ") if isinstance(x, str) else x
            )
            motif_dict = {"headers": motif_columns,
                               "rows": motif_consequences_dt.to_dict(orient="records")}

            # Get closest genes
            info_pd["DISTANCE"] = info_pd["DISTANCE"].replace('', pd.NA)
            float_distances = info_pd.dropna(subset = ["DISTANCE"])
            float_distances["DISTANCE"] = float_distances["DISTANCE"].astype("float")
            #float_distances_protein_coding = float_distances.loc[float_distances["BIOTYPE"] == "protein_coding"]
            if not float_distances.empty:
                closest_gene = [float_distances.loc[float_distances["DISTANCE"].astype("float").idxmin(), "SYMBOL"]]
            else:
                # take all genes (because no distance given)
                closest_gene = info_pd["SYMBOL"].drop_duplicates().tolist()
            rows+= info_dict # TODO: check if this even happens that you get two rows of a VCF file for a single position

        msc, msc_impact = get_most_severe(consequences)
        data = {
        "header": anno_columns,
        "rows": rows,
        "location": location,
        "ref": ref,
        "alt": alt,
        "external_ids": external_ids,
        "allele_frequencies": allele_frequencies,
        #"AF": AF,
        "msc": msc,
        "msc_impact": msc_impact,
        "transcript_consequences": transcript_dict if len(transcript_consequences_dt) > 0 else None,
        "regulatory_consequences": regulatory_dict if len(regulatory_consequences_dt) > 0 else None,
        "motif_consequences": motif_dict if len(motif_consequences_dt) > 0 else None,
        "closest_gene": closest_gene if len(closest_gene) > 0 else None,
        }
        return data
    except Exception as e:
        print("Error fetching data:" + str(e))
    return None  # no match found
