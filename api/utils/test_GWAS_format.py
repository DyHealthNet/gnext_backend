import gzip
import csv

file_path = "/nfs/data/Pan_UKBB/data/GWAS_stats_norm/biomarkers-30600-both_sexes-irnt_norm.tsv.bgz.gz"
with gzip.open(file_path, "rt") as f:
  reader = csv.DictReader(f, delimiter="\t")
  for row in reader:
      break
