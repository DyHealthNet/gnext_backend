#!/bin/bash

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /absolute/path/to/vep_dir GENOME_BUILD"
    exit 1
fi


VEP_DIR="$1"
GENOME_BUILD="$2"
VEP_CACHE_DIR="$VEP_DIR/vep_data"

TARGET_DIR=$(find "$VEP_CACHE_DIR/homo_sapiens" -maxdepth 1 -type d -name "*${GENOME_BUILD}" 2>/dev/null || true)

if [ -n "$TARGET_DIR" ]; then
    echo "VEP cache already exists for genome build '$GENOME_BUILD' in: $TARGET_DIR"
    echo "Skipping VEP setup steps."
    exit 0
fi

echo "No existing VEP cache found for $GENOME_BUILD. Proceeding with setup..."

echo "Creating VEP cache directory: $VEP_CACHE_DIR"
mkdir -p "$VEP_CACHE_DIR"

echo "Pulling VEP Docker image..."
docker pull ensemblorg/ensembl-vep

echo "Fix permissions"
chmod -R 777 "$VEP_CACHE_DIR"

echo "Running INSTALL.pl to download VEP cache..."
docker run -t -i \
  -v "$VEP_CACHE_DIR":/opt/vep/.vep \
  -e HOME=/opt/vep \
  --workdir /opt/vep \
  --user $(id -u):$(id -g) \
  ensemblorg/ensembl-vep \
  INSTALL.pl -a cf -s homo_sapiens -y "$GENOME_BUILD"

echo "Done! VEP cache is now prepared in: $VEP_CACHE_DIR"