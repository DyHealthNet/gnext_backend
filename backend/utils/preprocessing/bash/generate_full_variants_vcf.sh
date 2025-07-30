#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# Usage: ./generate_full_variants_vcf.sh /path/to/GWAS_dirs /path/to/output.vcf NUM_JOBS

# Print time
echo "Starting VCF generation at $(date)"

DATA_DIR="$1"
OUTPUT_FILE="$2"
NUM_JOBS="$3"

if [[ ! -d "$DATA_DIR" ]]; then
  echo "Usage: $0 /path/to/tsv_directory /path/to/output.vcf [NUM_JOBS]" >&2
  exit 1
fi

OUT_DIR=$(dirname "$OUTPUT_FILE")
[[ -d "$OUT_DIR" ]] || { echo "Error: $OUT_DIR does not exist." >&2; exit 1; }

# write header
{
  printf '##fileformat=VCFv4.2\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
} > "$OUTPUT_FILE"

# Check if pigz is available
if command -v pigz &>/dev/null; then
  USE_PIGZ=true
else
  USE_PIGZ=false
fi

# process everything in parallel, stream into a single sort/uniq
find "$DATA_DIR" -maxdepth 1 -type f -name '*.gz' \
  | parallel -j"$NUM_JOBS" --no-notice '
      PHENO={/.}
      # use pigz if available, fall back to gzip
      if [[ "$USE_PIGZ" == true ]]; then
        pigz -dc "{}"
      else
        gzip -cd "{}"
      fi \
      | awk -v pheno="$PHENO" -F"\t" -v OFS="\t" "
          NR==1 {
            for (i=1; i<=NF; i++) {
              h=tolower(\$i)
              if (h==\"#chrom\"||h==\"chrom\") c=i
              if (h==\"pos\")                  p=i
              if (h==\"ref\")                  r=i
              if (h==\"alt\")                  a=i
            }
            if (!c||!p||!r||!a) {
              print \"ERROR: header missing columns in \" pheno > \"/dev/stderr\"
              exit 1
            }
            next
          }
          \$c&&\$p&&\$r&&\$a {
            print \$c, \$p, \".\", \$r, \$a, \".\", \".\", \".\"
          }
      "
' \
  | sort -k1,1 -k2,2n -S10G \
  | uniq \
  >> "$OUTPUT_FILE"

chmod 777 "$OUTPUT_FILE"
echo "Done. VCF written to: $OUTPUT_FILE"

# Print end time
echo "Finished VCF generation at $(date)"