import os

# Пути
PDB_DIR = "/home/gluck/mnt/pdb/"         # корень архива PDB
DATA_FILE = "data/data_SORTED.csv"                  # входной файл
RESULTS_DIR = "/home/gluck/enrich/data/results/"                       # выходная папка для JSON
TEMP_DIR = "/home/gluck/ramdisk/dbenrich26/"                       # временные PDB-фрагменты
UNIPROT_CACHE_DIR = "/home/gluck/meta/uniprot/"     # кэш UniProt

# Управление обработкой
START_INDEX = 0            # индекс строки, с которой начинать
STOP_INDEX = 10000          # индекс строки, на которой останавливаться (None = до конца)
ONLY_NOT_PROCESSED = True  # пропускать уже обработанные (по наличию JSON)
STOP_ON_ERROR = False      # True = падаем при ошибке; False = логируем и продолжаем

# Параметры обработки
NUM_WORKERS = None         # None = os.cpu_count()
MISSING_RESIDUE_THRESHOLD = 0.10  # лимит пропущенных residue (10%)
COIL_TAIL_MAX_LEN = 3      # coil на конце <= этой длины считается "хвостом" и не учитывается

# Водородные связи: добавление атомов водорода перед расчётом
ADD_HYDROGENS = True       # True = добавлять H через obabel перед calcHydrogenBonds
OBABEL_BIN = "obabel"     # путь к бинарнику obabel (или просто "obabel" если в PATH)

# Параметры агрегатора
AGGREGATOR_START_INDEX = 0
AGGREGATOR_STOP_INDEX = None

# Кириллица допускается в комментариях
