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
load_env() {
    local env_file="../.env"
    if [ -f "$env_file" ]; then
        echo "Loading environment from $env_file"
        # Read only simple key=value pairs, ignore complex multi-line values
        while IFS='=' read -r key value; do
            # Skip empty lines, comments, and malformed lines
            [[ -z "$key" || "$key" =~ ^[[:space:]]*# || "$key" =~ [[:space:]] ]] && continue
            # Only process valid environment variable names
            if [[ "$key" =~ ^[A-Z_][A-Z0-9_]*$ ]]; then
                # Remove quotes if present and export
                value="${value%\"}"
                value="${value#\"}"
                export "$key=$value"
            fi
        done < <(grep -E '^[A-Z_][A-Z0-9_]*=' "$env_file" | head -20)
    else
        echo "Error: .env file not found at $env_file"
        exit 1
    fi
}

load_env

# --- Configuration ---
# Use defaults for development (VITE_TYPESENSE_HOST defaults to localhost, not from .env)
TYPESENSE_HOST="${VITE_TYPESENSE_HOST:-localhost}"
CONTAINER_NAME="gnext-typesense-${STUDY_ACRONYM:-dev}"
PORT="${VITE_TYPESENSE_PORT:-8108}"
API_KEY="${VITE_TYPESENSE_KEY}"
VOLUME_DIR="${TYPESENSE_VOLUME_DIR:-$(realpath ../typesense_volume)}"
IMAGE="typesense/typesense:29.0.rc15"

echo "================================"
echo "Typesense Development Setup"
echo "================================"
echo "Container: ${CONTAINER_NAME}"
echo "Host: ${TYPESENSE_HOST}"
echo "Port: ${PORT}"
echo "Volume: ${VOLUME_DIR}"
echo "Image: ${IMAGE}"
echo "Study: ${STUDY_ACRONYM:-dev}"
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
