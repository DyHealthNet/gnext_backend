#!/bin/bash

# Usage: ./run_vep.sh input.vcf output.vcf.bgz /absolute/path/to/vep_data genome_build window_up window_down num_jobs

set -euo pipefail

INPUT_VCF="$1"
OUTPUT_VCF="$2"
VEP_DIR="$3"
GENOME_BUILD=$4
WINDOW_UP=$5
WINDOW_DOWN=$6
NUM_JOBS="${7:-$(nproc)}"

VEP_CACHE_DIR="$VEP_DIR/vep_data/"

echo "Starting VEP run at $(date)"

echo $VEP_CACHE_DIR

docker run \
  -v "$VEP_DIR":/data \
  ensemblorg/ensembl-vep \
  vep \
    --input_file /data/"$INPUT_VCF" \
    --output_file /data/"$OUTPUT_VCF" \
    --cache \
    --offline \
    --dir_cache /data/vep_data \
    --compress_output bgzip \
    --fork "$NUM_JOBS" \
    --everything \
    --vcf \
    --species homo_sapiens \
    --assembly "$GENOME_BUILD" \
    --distance "$WINDOW_UP","$WINDOW_DOWN"

echo "VEP completed at $(date)"

tabix -p vcf "$VEP_DIR"/"$OUTPUT_VCF"