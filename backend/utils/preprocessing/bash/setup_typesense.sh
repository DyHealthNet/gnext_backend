#!/bin/bash

set -euo pipefail

# Check if a container named "typesense-server" is running
if docker ps --filter "name=dyhealthnetlight-typesense" --filter "status=running" | grep dyhealthnetlight-typesense > /dev/null; then
    echo "Typesense-server container is already running."
else
    echo "Typesense-server container is not running. Starting setup..."

    echo "Pulling Typesense image..."
    docker pull typesense/typesense:29.0.rc15

    echo "Creating Docker volume GWAS_typesense (if not exists)..."
    docker volume create GWAS_typesense

    echo "Starting typesense-server container..."
    docker run -d \
        --name dyhealthnetlight-typesense \
        -p 8108:8108 \
        -v GWAS_typesense:/data \
        typesense/typesense:29.0.rc15 \
        --data-dir /data \
        --backend \
        --api-key=xyz \
        --enable-cors

    echo "Typesense-server container started."
fi