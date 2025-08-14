# syntax=docker/dockerfile:1
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# System deps (add others you actually need, e.g. libpq-dev already here)
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gcc curl netcat-traditional libpq-dev \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Faster layer caching for deps
COPY requirements.txt .
RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# App source
COPY . .

# Make entrypoint executable
RUN chmod 0755 /app/docker-entrypoint.sh

# IMPORTANT: choose a user that can read/write your bind-mounts
# If your host user is UID 1000, this avoids bind-mount permission issues
RUN useradd -u 1000 -m django || true
USER 1000:1000

EXPOSE 8000
ENV DJANGO_SETTINGS_MODULE=dyhealthnetlight.settings

ENTRYPOINT ["/app/docker-entrypoint.sh"]
CMD ["gunicorn", "dyhealthnetlight.wsgi:application", "--bind", "0.0.0.0:8000", "--workers", "3"]
