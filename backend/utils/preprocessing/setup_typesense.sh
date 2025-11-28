#!/bin/bash
set -euo pipefail

# Usage:
#   ./setup_typesense.sh <container_name> <port> <api_key>
# Example:
#   ./setup_typesense.sh dyhealthnetlight-typesense-studyA 8109 xyz123

# --- Arguments ---
CONTAINER_NAME=$1
PORT=$2
API_KEY=$3
VOLUME_DIR=$4

# --- Static config ---
IMAGE="typesense/typesense:29.0.rc15"
VOLUME_NAME="${CONTAINER_NAME}_data"


# Check if the container is already running
if docker ps --filter "name=${CONTAINER_NAME}" --filter "status=running" | grep -q "${CONTAINER_NAME}"; then
    echo "Typesense container '${CONTAINER_NAME}' is already running."
    exit 0
fi

echo "Starting setup for Typesense container '${CONTAINER_NAME}'..."

# Pull latest image
echo "Pulling image '${IMAGE}'..."
docker pull "${IMAGE}"

# Create directory for volume
echo "Creating data directory at '${VOLUME_DIR}'..."
mkdir -p "${VOLUME_DIR}"

# Create volume if not exists
echo "Creating Docker volume '${VOLUME_NAME}' at '${VOLUME_DIR}'..."
docker volume create \
    --driver local \
    --opt type=none \
    --opt device="${VOLUME_DIR}" \
    --opt o=bind \
    "${VOLUME_NAME}" > /dev/null

# Run container
echo "Launching Typesense container '${CONTAINER_NAME}'..."
docker run -d \
    --name "${CONTAINER_NAME}" \
    -p "${PORT}:8108" \
    -v "${VOLUME_NAME}":/data \
    "${IMAGE}" \
    --data-dir /data \
    --api-key="${API_KEY}" \
    --enable-cors

echo "✅ Typesense container '${CONTAINER_NAME}' started successfully on port ${PORT}."
echo "   Volume: ${VOLUME_NAME}"