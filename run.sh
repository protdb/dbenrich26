#!/usr/bin/env bash
set -euo pipefail

IMAGE="${IMAGE:-glucksistemi/psskb_enricher}"
TAG="${TAG:-latest}"

CMD="${1:-processor.py}"

docker run --rm \
  -e PDB_DIR=/pdb \
  -e DATA_FILE=/data/data_SORTED.csv \
  -e RESULTS_DIR=/results \
  -e TEMP_DIR=/tmp/dbenrich26 \
  -e UNIPROT_CACHE_DIR=/home/gluck/meta/uniprot \
  -e STRIDE_BIN=/usr/local/bin/stride \
  -e START_INDEX="${START_INDEX:-1000}" \
  -e STOP_INDEX="${STOP_INDEX:-2000}" \
  -e NUM_WORKERS="${NUM_WORKERS:-}" \
  -e ONLY_NOT_PROCESSED="${ONLY_NOT_PROCESSED:-true}" \
  -e STOP_ON_ERROR="${STOP_ON_ERROR:-false}" \
  -e ADD_HYDROGENS="${ADD_HYDROGENS:-true}" \
  -e OBABEL_BIN="${OBABEL_BIN:-obabel}" \
  -e AGGREGATOR_START_INDEX="${AGGREGATOR_START_INDEX:-0}" \
  -e AGGREGATOR_STOP_INDEX="${AGGREGATOR_STOP_INDEX:-}" \
  -v "${PDB_VOLUME:-/home/gluck/mnt/pdb}:/pdb" \
  -v "${DATA_VOLUME:-$(pwd)/data}:/data" \
  -v "${RESULTS_VOLUME:-/home/gluck/enrich/data/results}:/results" \
  -v "${UNIPROT_VOLUME:-/home/gluck/meta/uniprot}:/uniprot_cache" \
  "${IMAGE}:${TAG}" \
  "$CMD"
