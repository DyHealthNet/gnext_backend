#!/usr/bin/env sh
set -e

# Optional: wait for Typesense on the Docker network
if [ -n "${VITE_TYPESENSE_HOST}" ] && [ -n "${VITE_TYPESENSE_PORT}" ]; then
  echo "Waiting for Typesense at ${VITE_TYPESENSE_HOST}:${VITE_TYPESENSE_PORT}..."
  for i in $(seq 1 60); do
    if nc -z "${VITE_TYPESENSE_HOST}" "${VITE_TYPESENSE_PORT}"; then
      echo "Typesense is up."
      break
    fi
    sleep 1
  done
fi

# Django housekeeping
python manage.py collectstatic --noinput || true
python manage.py migrate --noinput

# Hand off to CMD
exec "$@"
