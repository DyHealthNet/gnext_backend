#!/bin/bash

# Usage: ./run_vep.sh /absolute/path/to/input.vcf /absolute/path/to/output.vcf.bgz

set -euo pipefail

INPUT_VCF="$1"
OUTPUT_VCF="$2"
GENOME_BUILD=$3

if [[ -z "$INPUT_VCF" || ! -f "$INPUT_VCF" ]]; then
  echo "Error: Input VCF '$INPUT_VCF' does not exist."
  exit 1
fi

# Check parent directory of output file exists
OUT_DIR=$(dirname "$OUTPUT_FILE")
if [[ ! -d "$OUT_DIR" ]]; then
  echo "Error: Output directory '$OUT_DIR' does not exist."
  exit 1
fi

echo "Starting VEP run at $(date)"

docker run \
  -v /nfs/data/Pan_UKBB/data:/data \
  ensemblorg/ensembl-vep \
  vep \
    --input_file "$INPUT_VCF" \
    --output_file "$OUTPUT_VCF" \
    --cache \
    --offline \
    --compress_output bgzip \
    --fork 4 \
    --everything \
    --vcf \
    --species homo_sapiens \
    --assembly "$GENOME_BUILD"

echo "VEP completed at $(date)"

tabix -p vcf "$OUTPUT_VCF"