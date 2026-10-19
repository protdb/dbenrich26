import os


def _env_int(key: str, default: int) -> int:
    val = os.environ.get(key, '')
    return int(val) if val else default


def _env_float(key: str, default: float) -> float:
    val = os.environ.get(key, '')
    return float(val) if val else default


def _env_bool(key: str, default: bool) -> bool:
    val = os.environ.get(key, '').lower()
    if val in ('1', 'true', 'yes'):
        return True
    if val in ('0', 'false', 'no'):
        return False
    return default


def _env_optional_int(key: str) -> int | None:
    val = os.environ.get(key, '')
    return int(val) if val else None


PDB_DIR            = os.environ.get('PDB_DIR',            '/home/gluck/mnt/pdb/')
DATA_FILE           = os.environ.get('DATA_FILE',           'data/data_SORTED.csv')
RESULTS_DIR         = os.environ.get('RESULTS_DIR',         '/home/gluck/enrich/ext_result/results/')
TEMP_DIR            = os.environ.get('TEMP_DIR',            '/home/gluck/ramdisk/dbenrich26/')
UNIPROT_CACHE_DIR   = os.environ.get('UNIPROT_CACHE_DIR',   '/home/gluck/meta/uniprot/')

START_INDEX         = _env_int('START_INDEX',          0)
STOP_INDEX          = _env_optional_int('STOP_INDEX')
ONLY_NOT_PROCESSED  = _env_bool('ONLY_NOT_PROCESSED',  True)
STOP_ON_ERROR       = _env_bool('STOP_ON_ERROR',       False)

NUM_WORKERS             = _env_optional_int('NUM_WORKERS')
MISSING_RESIDUE_THRESHOLD = _env_float('MISSING_RESIDUE_THRESHOLD', 0.10)
COIL_TAIL_MAX_LEN        = _env_int('COIL_TAIL_MAX_LEN', 3)

ADD_HYDROGENS  = _env_bool('ADD_HYDROGENS', True)
OBABEL_BIN     = os.environ.get('OBABEL_BIN', 'obabel')

STRIDE_BIN     = os.environ.get('STRIDE_BIN', './stride')

AGGREGATOR_START_INDEX = _env_int('AGGREGATOR_START_INDEX', 1000000)
AGGREGATOR_STOP_INDEX  = _env_int('AGGREGATOR_STOP_INDEX', 1300000)
