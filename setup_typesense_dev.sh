#!/bin/bash
set -euo pipefail

# Development script to setup Typesense container
# Usage: ./setup_typesense_dev.sh [--pull]

# --- Parse arguments ---
PULL_IMAGE=false
if [[ "${1:-}" == "--pull" ]]; then
    PULL_IMAGE=true
fi

# Load environment variables from .env file
if [ -f ../.env ]; then
    export $(grep -v '^#' ../.env | xargs)
else
    echo "Error: .env file not found in parent directory"
    exit 1
fi

# --- Configuration ---
CONTAINER_NAME="${VITE_TYPESENSE_HOST}"
PORT="${VITE_TYPESENSE_PORT:-8108}"
API_KEY="${VITE_TYPESENSE_KEY}"
VOLUME_DIR="${TYPESENSE_DATA_DIR:-../typesense_volume}"
IMAGE="typesense/typesense:29.0.rc15"

echo "================================"
echo "Typesense Development Setup"
echo "================================"
echo "Container: ${CONTAINER_NAME}"
echo "Port: ${PORT}"
echo "Volume: ${VOLUME_DIR}"
echo "Image: ${IMAGE}"
echo "================================"

# Check if the container is already running
if docker ps --filter "name=${CONTAINER_NAME}" --filter "status=running" | grep -q "${CONTAINER_NAME}"; then
    echo "✅ Typesense container '${CONTAINER_NAME}' is already running."
    exit 0
fi

# Pull image if requested or if it doesn't exist
if [ "$PULL_IMAGE" = true ]; then
    echo "Pulling image '${IMAGE}'..."
    docker pull "${IMAGE}"
elif ! docker image inspect "${IMAGE}" > /dev/null 2>&1; then
    echo "Image '${IMAGE}' not found locally. Pulling..."
    docker pull "${IMAGE}"
else
    echo "Using existing local image '${IMAGE}'"
fi

# Create directory for volume
echo "Creating data directory at '${VOLUME_DIR}'..."
mkdir -p "${VOLUME_DIR}"

# Check if container exists but is stopped
if docker ps -a --filter "name=${CONTAINER_NAME}" | grep -q "${CONTAINER_NAME}"; then
    echo "Container '${CONTAINER_NAME}' exists but is stopped. Starting it..."
    docker start "${CONTAINER_NAME}"
else
    # Run new container
    echo "Launching Typesense container '${CONTAINER_NAME}'..."
    docker run -d \
        --name "${CONTAINER_NAME}" \
        -p "${PORT}:8108" \
        -v "${VOLUME_DIR}":/data \
        "${IMAGE}" \
        --data-dir /data \
        --api-key="${API_KEY}" \
        --enable-cors
fi

echo "✅ Typesense container '${CONTAINER_NAME}' is running on port ${PORT}."
echo ""
echo "To initialize data, run:"
echo "  python init_typesense.py"
