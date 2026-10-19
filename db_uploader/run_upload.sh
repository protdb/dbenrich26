#!/usr/bin/env bash
set -euo pipefail

IMAGE="${IMAGE:-glucksistemi/psskb_uploader}"
TAG="${TAG:-latest}"

PG_DSN=postgresql://data_api:123QWErty@172.17.0.2:5432/data_api

DATA_VOLUME=/home/gluck/PycharmProjects/dbenrich26/data

if [[ -z "${PG_DSN:-}" ]]; then
  echo "ERROR: Переменная PG_DSN не задана"
  echo "Использование: PG_DSN=\"postgresql://user:pass@host:5432/dbname\" DATA_VOLUME=\"/path/to/csv\" ./run_upload.sh"
  exit 1
fi

docker build -t "${IMAGE}:${TAG}" .

docker run --rm \
  -e PG_DSN="${PG_DSN}" \
  -e BATCH_SIZE="${BATCH_SIZE:-500}" \
  -e WRITE_MODE="${WRITE_MODE:-update}" \
  -v "${DATA_VOLUME:-$(pwd)/data}:/data" \
  "${IMAGE}:${TAG}"
