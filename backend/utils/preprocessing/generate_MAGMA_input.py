from decouple import config
import gzip
from collections import defaultdict
import time

without_gene_terms = ["regulatory_region_variant", "TF_binding_site_variant", "intergenic_variant", "intron_variant", "upstream_gene_variant", "downstream_gene_variant"]

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

def extract_protein_coding_rsids(vcf_path):
    csq_fields = extract_csq_fields(vcf_path)
    gene_to_rsids = defaultdict(set)
    i = 1
    with gzip.open(vcf_path, 'rt') as f:
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
            i+=1
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
                    gene_to_rsids[gene].add(rsid)

                elif any(c in without_gene_terms for c in consequences):
                    try:
                        distance = int(csq_dict.get("DISTANCE", "999999"))  # fallback to large if missing
                    except ValueError:
                        continue
                    if distance < window_size:
                        gene_to_rsids[gene].add(rsid)
            if i % 100000 == 0:
                print(f"Processed {i} lines...")
                print(f"Current gene to rsid mapping size: {len(gene_to_rsids)}")
    return gene_to_rsids


vep_anno_file = "/Users/lisiarend/Desktop/DyHealthNet/DyHealthNetLight/data/annotated_variants.vcf.bgz"
window_size = 1000  # Size of the window for tabix queries

# Stop the time
start_time = time.time()
mapping = extract_protein_coding_rsids(vep_anno_file)
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Time taken to process the VCF file: {elapsed_time:.2f} seconds")
print(f"Total unique genes found: {len(mapping)}")
