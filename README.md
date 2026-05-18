# PSSKB Enricher

Pipeline for enriching PDB motif data with SASA, secondary structure, inter-residue
interactions, and UniProt metadata.

## Local development

```bash
pip install poetry==2.2.1
poetry install
```

Run the processor:

```bash
poetry run python processor.py
```

Run the aggregator (after processing):

```bash
poetry run python aggregator.py
```

## Docker

```bash
docker pull glucksistemi/psskb_enricher
```

### Processor

```bash
docker run --rm \
  -e START_INDEX=0 \
  -e STOP_INDEX=1000 \
  -e ONLY_NOT_PROCESSED=true \
  -v /path/to/pdb:/pdb \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  -v /path/to/uniprot_cache:/uniprot_cache \
  glucksistemi/psskb_enricher processor.py
```

### Aggregator

```bash
docker run --rm \
  -e AGGREGATOR_START_INDEX=0 \
  -e AGGREGATOR_STOP_INDEX=1000 \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  glucksistemi/psskb_enricher aggregator.py
```

### Build & push

```bash
./build_and_push.sh          # latest
./build_and_push.sh v1.0.0   # specific tag
```

## Environment variables

### Paths

| Variable | Default | Description |
|---|---|---|
| `PDB_DIR` | `/home/gluck/mnt/pdb/` | PDB archive root |
| `DATA_FILE` | `data/data_SORTED.csv` | Input CSV path |
| `RESULTS_DIR` | `/home/gluck/enrich/data/results/` | JSON results output |
| `TEMP_DIR` | `/home/gluck/ramdisk/dbenrich26/` | Temporary PDB fragments |
| `UNIPROT_CACHE_DIR` | `/home/gluck/meta/uniprot/` | UniProt metadata cache |
| `STRIDE_BIN` | `./stride` | Path to STRIDE binary |

### Processing control

| Variable | Default | Description |
|---|---|---|
| `START_INDEX` | `0` | First CSV row to process |
| `STOP_INDEX` | *(end of file)* | Last CSV row + 1 |
| `NUM_WORKERS` | *(CPU count)* | Number of parallel workers |
| `ONLY_NOT_PROCESSED` | `true` | Skip already-processed entries |
| `STOP_ON_ERROR` | `false` | Crash on first error |

### Interactions

| Variable | Default | Description |
|---|---|---|
| `ADD_HYDROGENS` | `true` | Add H atoms via obabel before HB calc |
| `OBABEL_BIN` | `obabel` | Path to obabel binary |

### Aggregator

| Variable | Default | Description |
|---|---|---|
| `AGGREGATOR_START_INDEX` | `0` | First row to aggregate |
| `AGGREGATOR_STOP_INDEX` | *(end)* | Last row + 1 |
