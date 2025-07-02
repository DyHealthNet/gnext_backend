#!/bin/bash

# Usage: ./generate_full_variants_vcf.sh /path/to/tsv_directory /path/to/output.vcf

set -euo pipefail

DATA_DIR="$1"
OUTPUT_FILE="$2"

if [[ -z "$DATA_DIR" || ! -d "$DATA_DIR" ]]; then
  echo "Usage: $0 /path/to/GWAS_stats_norm /path/to/output.vcf"
  exit 1
fi

# Check parent directory of output file exists
OUT_DIR=$(dirname "$OUTPUT_FILE")
if [[ ! -d "$OUT_DIR" ]]; then
  echo "Error: Output directory '$OUT_DIR' does not exist."
  exit 1
fi

TMP_DIR="$OUT_DIR/tmp_variants"
RAW_BODY="$TMP_DIR/variants_raw_body.tsv"

echo "Creating working directory: $TMP_DIR"
mkdir -p "$TMP_DIR"
rm -f "$RAW_BODY" "$OUTPUT_FILE"

# Save VCF header
echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$OUTPUT_FILE"

echo "Processing .gz files in $DATA_DIR..."

for f in "$DATA_DIR"/*.gz; do
  [[ -e "$f" ]] || continue
  PHENO=$(basename "$f" .gz)

  echo "Extracting variants from $PHENO..."

  gzcat "$f" | awk -v pheno="$PHENO" -F'\t' '
    BEGIN {OFS="\t"}
    NR==1 {
      for (i=1; i<=NF; i++) {
        if (tolower($i)=="#chrom" || tolower($i)=="chrom") c=i;
        if (tolower($i)=="pos") p=i;
        if (tolower($i)=="ref") r=i;
        if (tolower($i)=="alt") a=i;
      }
      if (!c || !p || !r || !a) {
        print "ERROR: Required columns not found in header of " pheno > "/dev/stderr"
        exit 1
      }
      next
    }
    $c && $p && $r && $a {
      print $c, $p, ".", $r, $a, ".", ".", "."
    }
  ' >> "$RAW_BODY"
done

echo "Deduplicating and sorting variants..."
sort -k1,1 -k2,2n "$RAW_BODY" | uniq >> "$OUTPUT_FILE"

echo "Done."
echo "VCF written to: $OUTPUT_FILE"