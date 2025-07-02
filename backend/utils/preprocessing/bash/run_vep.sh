#!/bin/bash

# Usage: ./run_vep.sh input.vcf output.vcf.bgz /absolute/path/to/vep_data genome_build

set -euo pipefail

INPUT_VCF="$1"
OUTPUT_VCF="$2"
VEP_CACHE_DIR="$3"
GENOME_BUILD=$4

echo "Starting VEP run at $(date)"

docker run \
  -v "$VEP_CACHE_DIR":/data \
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