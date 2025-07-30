#!/bin/bash

# Usage: ./generate_full_variants_vcf.sh /path/to/tsv_directory /path/to/output.vcf NUM_JOBS

set -euo pipefail

DATA_DIR="$1"
OUTPUT_FILE="$2"
NUM_JOBS="$3"

if [[ -z "$DATA_DIR" || ! -d "$DATA_DIR" ]]; then
  echo "Usage: $0 /path/to/GWAS_stats_norm /path/to/output.vcf [NUM_JOBS]"
  exit 1
fi

OUT_DIR=$(dirname "$OUTPUT_FILE")
if [[ ! -d "$OUT_DIR" ]]; then
  echo "Error: Output directory '$OUT_DIR' does not exist."
  exit 1
fi

TMP_DIR="$OUT_DIR/tmp_variants"
mkdir -p "$TMP_DIR"
trap 'rm -rf "$TMP_DIR"' EXIT

echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$OUTPUT_FILE"

process_file() {
  FILE="$1"
  PHENO=$(basename "$FILE" .gz)
  TMP_OUT="$TMP_DIR/variants_${PHENO}.tsv"

  gzip -cd "$FILE" | awk -v pheno="$PHENO" -F'\t' '
    BEGIN { OFS = "\t" }
    NR==1 {
      for (i=1; i<=NF; i++) {
        if (tolower($i)=="#chrom"||tolower($i)=="chrom") c=i
        if (tolower($i)=="pos")                   p=i
        if (tolower($i)=="ref")                   r=i
        if (tolower($i)=="alt")                   a=i
      }
      if (!c || !p || !r || !a) {
        print "ERROR: Missing columns in header of " pheno > "/dev/stderr"
        exit 1
      }
      next
    }
    $c && $p && $r && $a {
      print $c, $p, ".", $r, $a, ".", ".", "."
    }
  ' > "$TMP_OUT"
    echo "Processed ${PHENO}" >&2
}

export -f process_file
export TMP_DIR

echo "Extracting variants using $NUM_JOBS parallel jobs..."
find "$DATA_DIR" -maxdepth 1 -type f -name '*.gz' \
  | parallel --line-buffer -j"$NUM_JOBS" process_file {}

echo "Sorting and deduplicating..."
cat "$TMP_DIR"/variants_*.tsv \
  | sort -k1,1 -k2,2n -S10G \
  | uniq >> "$OUTPUT_FILE"

# change permissions of the output file
chmod 777 "$OUTPUT_FILE"

echo "Done. VCF written to: $OUTPUT_FILE"