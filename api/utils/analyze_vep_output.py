import gzip
from decouple import config
import pysam
from converters import convert_variant_id

def extract_variant_annotation(variant_id):
    vep_anno_file = config('VEP_ANNO_FILE')

    # Extract GWAS results of variant from phenotype file via tabix
    tabix_file = pysam.TabixFile(vep_anno_file)
    # Get variant information

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
        rows = []
        for row in tabix_file.fetch(chr, pos - 1, pos):
            row = row.split("\t")
            row = dict(zip(columns, row))
            info = row["INFO"].replace("CSQ=", "").split(",")
            info_dict = [dict(zip(anno_columns, i.split("|"))) for i in info]
            rows = rows + info_dict
        data = {
            "header": anno_columns,
            "rows": rows,
        }
        return data
    except Exception as e:
        print("Error fetching data:" + str(e))
    return None  # no match found


# vep_anno_file = config('VEP_ANNO_FILE')
# n = 0
# with gzip.open(vep_anno_file, 'rt') as vcf:
#     for line in vcf:
#         if line.startswith("#CHROM"):
#             header = line.strip().split("\t")
#             header = [h.replace("#", "") for h in header]
#             print(header)
#             continue
#         if line.startswith("#"):
#             continue
#         fields = line.strip().split("\t")
#         variant_dict = dict(zip(header, fields))
#         info = variant_dict["INFO"].split(",")
#         rs_ids = list(set([i.split("|")[17] for i in info]))
#         rs_ids = ", ".join(rs_ids)
#         print(rs_ids)
#         n+= 1
#         if n == 10:
#             break


t = extract_variant_annotation("12_91999816_C/G")
print(t)