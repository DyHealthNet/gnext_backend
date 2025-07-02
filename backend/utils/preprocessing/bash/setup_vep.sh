#!/bin/bash

set -euo pipefail

VEP_CACHE_DIR="$1"
GENOME_BUILD="$2"

echo "Creating VEP cache directory: $VEP_CACHE_DIR"
mkdir -p "$VEP_CACHE_DIR"

echo "Pulling VEP Docker image..."
docker pull ensemblorg/ensembl-vep

echo "Running INSTALL.pl to download VEP cache..."
docker run -t -i \
  -v "$VEP_CACHE_DIR":/opt/vep/.vep \
  ensemblorg/ensembl-vep \
  INSTALL.pl -a cf -s homo_sapiens -y "$GENOME_BUILD"

echo "Done! VEP cache is now prepared in: $VEP_CACHE_DIR"
echo "You can now start your services using docker-compose up."