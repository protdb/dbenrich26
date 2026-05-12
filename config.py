import os

# Пути
PDB_DIR = "/home/gluck/mnt/Data2/PDB/pdb/"         # корень архива PDB
DATA_FILE = "data/data_SORTED.csv"                  # входной файл
RESULTS_DIR = "data/results/"                       # выходная папка для JSON
TEMP_DIR = "/tmp/dbenrich26/"                       # временные PDB-фрагменты
UNIPROT_CACHE_DIR = "/home/gluck/meta/uniprot/"     # кэш UniProt

# Управление обработкой
START_INDEX = 0            # индекс строки, с которой начинать
STOP_INDEX = None          # индекс строки, на которой останавливаться (None = до конца)
ONLY_NOT_PROCESSED = True  # пропускать уже обработанные (по наличию JSON)
STOP_ON_ERROR = False      # True = падаем при ошибке; False = логируем и продолжаем

# Параметры обработки
NUM_WORKERS = None         # None = os.cpu_count()
MISSING_RESIDUE_THRESHOLD = 0.10  # лимит пропущенных residue (10%)
COIL_TAIL_MAX_LEN = 3      # coil на конце <= этой длины считается "хвостом" и не учитывается

# Параметры агрегатора
AGGREGATOR_START_INDEX = 0
AGGREGATOR_STOP_INDEX = None

# Кириллица допускается в комментариях
