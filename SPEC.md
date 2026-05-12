# SPEC: dbenrich26 — Обогащение базы белковых структурных мотивов

## Обзор проекта

Набор скриптов для обработки 1.2M строк данных о белковых структурных мотивах из PDB.  
Входной файл: `data/data_SORTED.csv` (отсортирован по `experiment_id`, `chain_id`, `motif_start`).  
Выходной формат: CSV с теми же полями, что входной, плюс новое поле `gravy`. Пересчитываемые поля перезаписываются.  
Промежуточное хранилище: JSON-файлы по одному на мотив.

---

## Архитектура: файлы проекта

```
dbenrich26/
├── config.py               # Конфигурационные константы
├── stride                  # Бинарник STRIDE (уже есть)
├── processor.py            # Основной скрипт обработки (multiprocessing)
├── stride_tool.py          # Модуль работы со STRIDE
├── sasa_tool.py            # Модуль расчёта SASA (FreeSASA)
├── interactions_tool.py    # Модуль расчёта взаимодействий (ProDy)
├── structure_tool.py       # Модуль: загрузка PDB, вырезание мотива, Rg, масса, FASTA
├── uniprot_tool.py         # Модуль: получение данных UniProt
├── aggregator.py           # Скрипт конвертации JSON → CSV
└── data/
    ├── data_SORTED.csv
    └── results/            # Папка для JSON-результатов (задаётся в config)
        ├── 0000/
        ├── 1000/
        └── ...
```

---

## Выходные данные

### Поля в JSON / итоговом CSV

Только **пересчитанные** поля (+ новое `gravy`):

| Поле | Тип | Описание |
|---|---|---|
| `apfid` | str | Идентификатор (ключ) |
| `chain_len` | int | Длина цепи (только белковые residue из PDB) |
| `motif_len` | int | `motif_end - motif_start` |
| `fasta` | str | Последовательность мотива из PDB, нестандартные → `?` |
| `organism` | str | `up_meta.organism.scientific_name` |
| `organism_id` | str | `up_meta.organism.taxon_id` |
| `protein_name` | str | `up_meta.name` |
| `biological_process` | str | из UniProt keywords |
| `molecular_function` | str | из UniProt keywords |
| `rg` | float | Радиус гирации по Cα мотива |
| `sasa` | float | SASA мотива в изоляции (только его цепь) |
| `surface_k` | float | SASA мотива в контексте цепи / SASA в изоляции |
| `inner_hb` | int | H-связи внутри мотива (ProDy) |
| `inner_hp` | int | Гидрофобные взаимодействия внутри мотива (ProDy) |
| `inner_ds` | int | Дисульфидные мостики внутри мотива (ProDy) |
| `outer_hb` | int | H-связи мотив↔остаток цепи (ProDy) |
| `outer_hp` | int | Гидрофобные взаимодействия мотив↔остаток (ProDy) |
| `outer_ds` | int | Дисульфидные мостики мотив↔остаток (ProDy) |
| `motif_weight` | float | Молекулярная масса мотива (BioPython) |
| `protein_weight` | float | Масса из `up_meta.sequence` (BioPython) |
| `weight_percent` | float | `motif_weight / protein_weight * 100` |
| `protein_partners` | str | из UniProt |
| `has_binding_sites` | str | `"Yes"` / `"No"` |
| `binding_sites` | str | позиции сайтов связывания из UniProt |
| `helix_count` | int | Кол-во спиралей в мотиве (STRIDE, без хвостов) |
| `strand_count` | int | Кол-во тяжей (STRIDE, без хвостов) |
| `coil_count` | int | Кол-во coil-элементов (STRIDE, без хвостов) |
| `helix_length` | int | Суммарная длина спиралей |
| `strand_length` | int | Суммарная длина тяжей |
| `coil_length` | int | Суммарная длина coil-элементов |
| `gravy` | float | GRAVY-индекс мотива (шкала Kyte-Doolittle) |

### JSON для некорректных структур

```json
{
  "apfid": "...",
  "error": "описание причины"
}
```
Имя файла: `{apfid}__incorrect.json`

---

## Детальная декомпозиция задач

---

### ЗАДАЧА 0: Настройка зависимостей

**Файл:** `pyproject.toml`  
**Зависимостей нет**

Добавить в `pyproject.toml`:
- `biopython` (последняя стабильная)
- `freesasa`
- `prody`

Установить через `poetry add`.

---

### ЗАДАЧА 1: `config.py`

**Файл:** `config.py`  
**Зависит от:** —

Создать файл с константами:

```python
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
```

---

### ЗАДАЧА 2: `structure_tool.py`

**Файл:** `structure_tool.py`  
**Зависит от:** ЗАДАЧА 0 (biopython), ЗАДАЧА 1  
**Библиотеки:** `Bio.PDB`, `gzip`

Реализует:

#### 2.1 `get_pdb_path(experiment_id: str) -> str`
- Конструирует путь по шаблону: `{PDB_DIR}/{средние_2_символа}/pdb{experiment_id_lower}.ent.gz`
- Пример: `4E06` → `/home/gluck/mnt/Data2/PDB/pdb/e0/pdb4e06.ent.gz`

#### 2.2 `load_structure(experiment_id: str) -> Bio.PDB.Structure`
- Открывает `.ent.gz` через `gzip.open`
- Парсит через `Bio.PDB.PDBParser(QUIET=True)`
- Возвращает структуру (объект Structure)
- Бросает исключение если файл не найден

#### 2.3 `get_protein_residues(chain) -> list[Residue]`
- Фильтрует residue цепи: только стандартные и **модифицированные** аминокислоты
- Исключает: воду (`HOH`), лиганды, ионы (residue у которых `residue.id[0] != ' '` — hetflag)
- **Важно:** стандартный белковый остаток имеет `residue.id[0] == ' '`. Модифицированные также имеют `' '` в hetflag если они в пептидной цепи, иначе `'H_XXX'` — нужно проверять наличие атома `CA`.
  - Критерий включения: `residue.id[0] == ' '` ИЛИ (`residue.id[0].startswith('H_')` И наличие атома `'CA'`)

#### 2.4 `get_chain_len(chain) -> int`
- Возвращает `len(get_protein_residues(chain))`

#### 2.5 `extract_motif_residues(chain, motif_start: int, motif_end: int) -> list[Residue]`
- Из `get_protein_residues(chain)` берёт residue с `residue.id[1]` (sequence number) в диапазоне `[motif_start, motif_end]` включительно
- iCode игнорируем (`residue.id[2]`)
- Возвращает список Residue

#### 2.6 `validate_motif(motif_residues: list, motif_start: int, motif_end: int) -> None`
- Ожидаемая длина: `motif_end - motif_start`
- Реальная: `len(motif_residues)`
- Если `реальная < ожидаемая * (1 - MISSING_RESIDUE_THRESHOLD)` → бросает `MotifValidationError`
- Если список пуст → бросает `MotifValidationError`

#### 2.7 `save_motif_to_pdb(motif_residues: list, chain_id: str, out_path: str) -> None`
- Создаёт временный PDB-файл содержащий только residue мотива
- Использует `Bio.PDB.PDBIO` + `Select`-класс, фильтрующий по списку residue
- Путь: `{TEMP_DIR}/{apfid}.pdb`

#### 2.8 `save_chain_to_pdb(chain, out_path: str) -> None`
- Сохраняет всю одну цепь в PDB-файл (для расчёта SASA в контексте цепи)
- Путь: `{TEMP_DIR}/{apfid}_chain.pdb`

#### 2.9 `calc_rg(motif_residues: list) -> float`
- Радиус гирации по атомам Cα мотива
- Формула: `sqrt( sum((r_i - r_cm)^2) / N )`, где `r_cm` — центр масс Cα
- Бросает исключение если нет ни одного Cα

#### 2.10 `calc_motif_weight(motif_residues: list) -> float`
- Использует `Bio.SeqUtils.molecular_weight(seq, seq_type='protein')`
- Последовательность строим из residue (через `get_motif_fasta`)
- `X` (нестандартные) считаем как среднюю массу аминокислотного остатка ≈ 111.1 Да

#### 2.11 `get_motif_fasta(motif_residues: list) -> str`
- Для каждого residue пробуем `protein_letters_3to1[res.resname]`
- Если KeyError → `'?'`
- Возвращает строку

#### 2.12 `calc_protein_weight(sequence: str) -> float`
- Принимает строку FASTA полного белка из UniProt
- `X` → средняя масса остатка 111.1 Да
- Если `sequence` пустой или None → возвращает `0.0`

---

### ЗАДАЧА 3: `stride_tool.py`

**Файл:** `stride_tool.py`  
**Зависит от:** ЗАДАЧА 1, ЗАДАЧА 2 (пути к бинарнику STRIDE, путь к temp PDB)

Переписывает `examples/stride.py` с улучшениями.

#### 3.1 `run_stride(pdb_path: str) -> str | None`
- Запускает `./stride -f{out_path} {pdb_path}` через `subprocess`
- Выходной файл: `{pdb_path}.stride` (удаляется после парсинга)
- Возвращает raw-текст или None при ошибке
- После получения результата — удаляет `.stride`-файл

#### 3.2 `parse_stride_output(raw: str) -> list[dict]`
- Парсит ASG-строки STRIDE
- Поля: `amino`, `total_no` (sequence number), `chain`, `sec_letter`, `sec_name`
- Ремаппинг: `T`, `B`, `G` → `'C'` (coil)
- Возвращает список по всем residue

#### 3.3 `get_secondary_structure_elements(parsed: list[dict]) -> list[dict]`
- Схлопывает residue в элементы вторичной структуры
- Каждый элемент: `{'type': str, 'sec_letter': str, 'start': int, 'end': int, 'length': int}`
- `length = end - start + 1`

#### 3.4 `get_motif_ss_stats(elements: list[dict], motif_start: int, motif_end: int, coil_tail_max_len: int = COIL_TAIL_MAX_LEN) -> dict`
- Фильтрует элементы: только те, которые **полностью** входят в `[motif_start, motif_end]`
- Удаляет "хвостовые" coil-элементы:
  - Если **первый** элемент — `'C'` и `length <= coil_tail_max_len` → удаляем
  - Если **последний** элемент — `'C'` и `length <= coil_tail_max_len` → удаляем
  - (однократно для каждого конца, не рекурсивно)
- Считает и возвращает:
  ```python
  {
      'helix_count': int,   # элементы с sec_letter == 'H' или 'I'
      'strand_count': int,  # элементы с sec_letter == 'E'
      'coil_count': int,    # элементы с sec_letter == 'C'
      'helix_length': int,
      'strand_length': int,
      'coil_length': int,
  }
  ```

---

### ЗАДАЧА 4: `sasa_tool.py`

**Файл:** `sasa_tool.py`  
**Зависит от:** ЗАДАЧА 0 (freesasa), ЗАДАЧА 2 (структуры BioPython в памяти)

> **Важно:** FreeSASA принимает структуру напрямую из BioPython — не через путь к файлу.  
> Временный PDB-файл нужен **только** для STRIDE (который не умеет работать с объектами в памяти).  
> Для SASA используем объекты `Bio.PDB.Structure` / `Bio.PDB.Chain`, уже загруженные в память, чтобы избежать лишних I/O операций.

#### 4.1 `calc_sasa_isolated(motif_residues: list[Residue]) -> float`
- Принимает список residue мотива (BioPython объекты из памяти)
- Собирает временную структуру в памяти: создаёт `Bio.PDB.Structure` → `Model` → `Chain` с только этими residue
- Вычисляет SASA через `freesasa.calcBioPDB(tmp_structure)`
- Суммирует SASA по всем атомам полученной структуры
- Возвращает суммарное значение SASA в Å²

#### 4.2 `calc_sasa_in_chain(chain, motif_residue_ids: list[int]) -> float`
- Принимает объект `Bio.PDB.Chain` (вся цепь, уже в памяти) и список sequence numbers мотива
- Вычисляет SASA для полной цепи через `freesasa.calcBioPDB(chain_structure)`
- Суммирует SASA **только по атомам residue с `resseq` из `motif_residue_ids`**
- Возвращает суммарную SASA этих residue в контексте цепи

#### 4.3 `calc_surface_k(sasa_isolated: float, sasa_in_chain: float) -> float`
- `surface_k = sasa_in_chain / sasa_isolated`
- Если `sasa_isolated == 0` → возвращает `0.0`

---

### ЗАДАЧА 5: `interactions_tool.py`

**Файл:** `interactions_tool.py`  
**Зависит от:** ЗАДАЧА 0 (prody), ЗАДАЧА 2 (save_chain_to_pdb)

Использует ProDy для расчёта межостаточных взаимодействий.

#### 5.1 `calc_interactions(chain_pdb_path: str, motif_residue_ids: set[int]) -> dict`
- Загружает цепь через `prody.parsePDB(chain_pdb_path)`
- Вычисляет взаимодействия через `prody.Interactions()`:
  - `calcHydrogenBonds()` → hb
  - `calcHydrophobicContacts()` → hp  
  - `calcDisulfideBonds()` → ds
- Для каждого типа взаимодействия разбивает пары на inner/outer:
  - **inner**: оба residue (`resnum`) в `motif_residue_ids`
  - **outer**: один residue в `motif_residue_ids`, второй — нет
  - Пары с обоими вне мотива игнорируем
- Возвращает:
  ```python
  {
      'inner_hb': int, 'inner_hp': int, 'inner_ds': int,
      'outer_hb': int, 'outer_hp': int, 'outer_ds': int,
  }
  ```

---

### ЗАДАЧА 6: `uniprot_tool.py`

**Файл:** `uniprot_tool.py`  
**Зависит от:** ЗАДАЧА 1 (UNIPROT_CACHE_DIR)  
**Библиотеки:** `uniprot_meta_tool`

#### 6.1 `init_uniprot(cache_dir: str) -> None`
- Вызывает `set_config(meta_cache_dir=cache_dir)`

#### 6.2 `get_uniprot_data(uniprot_id: str) -> dict`
- Создаёт `UniprotData(uniprot_id)`
- Обрабатывает случаи: пустой `uniprot_id`, `ObsoleteException`, сетевые ошибки → возвращает `None`
- Извлекает и возвращает словарь:
  ```python
  {
      'protein_name': up.name,
      'biological_process': '; '.join(get_biological_processes(up)),
      'molecular_function': '; '.join(get_molecular_functions(up)),
      'protein_weight': up.sequence (передаётся в structure_tool.calc_protein_weight),
      'organism': up.organism.scientific_name if up.organism else None,
      'organism_id': up.organism.taxon_id if up.organism else None,
      'protein_partners': ...,    # см. логику ниже
      'has_binding_sites': ...,
      'binding_sites': ...,
  }
  ```

#### 6.3 Логика `protein_partners`

Готовый рецепт из `examples/get_up_additionals.py` (строки 68–77), адаптированный для нашего контекста.  
В качестве источника списка цепей структуры и их UniProt ID используем PDB API (`get_pdb_api_data.get_single`), как в примере:

```python
# plain_data — результат jmespath-запроса DATA_PLAIN к PDB API для данного experiment_id
# up_id — uniprot_id текущей записи (из CSV)
partner_uniprots = list(set([x['uniprot_id'] for x in plain_data if x['uniprot_id'] != up_id]))
if len(partner_uniprots) > 0:
    try:
        partner_ups = [UniprotData(x) for x in partner_uniprots if type(x) is str]
        protein_partners = [f"{x.gene} ({x.uniprot_id})" if x.gene else f"{x.uniprot_id}" for x in partner_ups]
    except TypeError:
        protein_partners = []
else:
    protein_partners = []
protein_partners = '; '.join(protein_partners)
```

- `plain_data` получается через `get_pdb_api_data.get_single(experiment_id)` + jmespath с запросом `DATA_PLAIN` (строка 33 примера)
- Ошибки запроса к PDB API → `protein_partners = ''`

#### 6.4 Логика `has_binding_sites` / `binding_sites`
- `features = up.get_position_features(motif_start, motif_end)`
- `binding_site_features = [x for x in features if x.feature_type == 'Binding site']`
- `has_binding_sites = "Yes" if len(...) > 0 else "No"`
- `binding_sites = '; '.join([str(x.start) if x.start == x.end else f'{x.start}-{x.end}' ...])`

#### 6.5 Fallback при недоступном UniProt
Если `uniprot_id` пустой, неизвестный, или запрос завершился ошибкой:
- `protein_name`, `biological_process`, `molecular_function`, `organism`, `organism_id` → оставляем из исходных данных CSV (не перезаписываем)
- `protein_weight = 0.0`, `weight_percent = 0.0`
- `protein_partners = ''`, `has_binding_sites = 'No'`, `binding_sites = ''`

---

### ЗАДАЧА 7: GRAVY-индекс

**Место:** вынести в `structure_tool.py` как функцию (или отдельный мини-модуль)  
**Зависит от:** ЗАДАЧА 2.11 (`get_motif_fasta`), константы из CornerPhobia

#### 7.1 `calc_gravy(fasta: str) -> float`
- Шкала Kyte-Doolittle (`KD_SCALE` из `CornerPhobia/corner_hydropathy/constants.py` — скопировать константы в проект)
- GRAVY = среднее значение гидрофобности по всем residue последовательности
- `?` (нестандартные) — пропускаем (не включаем в среднее)
- Если нет ни одного известного residue → возвращает `0.0`

---

### ЗАДАЧА 8: `processor.py` — основной скрипт

**Файл:** `processor.py`  
**Зависит от:** ЗАДАЧА 1, 2, 3, 4, 5, 6, 7

Основная логика обработки.

#### 8.1 Чтение входных данных

Файл `data_SORTED.csv` огромный (1.2M строк) — **не читаем его через pandas целиком**.  
Используем модуль `csv` и потоковый генератор групп.

**Генератор `iter_experiment_groups(filepath, start_index, stop_index)`:**

```python
import csv

def iter_experiment_groups(filepath, start_index=0, stop_index=None):
    """
    Читает CSV построчно, группируя строки по experiment_id.
    Файл должен быть отсортирован по experiment_id.
    Выдаёт (experiment_id, rows_list, global_indices_list) — по одной группе за раз.
    """
    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        current_group_id = None
        current_rows = []
        current_indices = []

        for abs_index, row in enumerate(reader):
            # Применяем срез START_INDEX..STOP_INDEX
            if abs_index < start_index:
                continue
            if stop_index is not None and abs_index >= stop_index:
                break

            exp_id = row['experiment_id']
            if exp_id != current_group_id:
                if current_group_id is not None:
                    yield (current_group_id, current_rows, current_indices)
                current_group_id = exp_id
                current_rows = []
                current_indices = []

            current_rows.append(row)
            current_indices.append(abs_index)

        if current_group_id is not None:
            yield (current_group_id, current_rows, current_indices)
```

- Каждый `yield` — одна задача для пула: `(experiment_id, rows_list, global_indices)`
- Генератор **не хранит в памяти весь файл** — только текущую группу
- Передаётся в `pool.imap()` (не `pool.map()`) — это позволяет не материализовывать все задачи сразу

#### 8.2 Структура папок результатов

- `RESULTS_DIR` создаётся если не существует
- Подпапки: `0000/`, `1000/`, `2000/` и т.д.
- Номер подпапки: `(global_row_index // 1000) * 1000` в нулях с форматом `{:04d}`
- `global_row_index` = позиция строки в срезе `START_INDEX..STOP_INDEX`

#### 8.3 `get_result_path(row_index: int, apfid: str, incorrect: bool = False) -> str`
- Подпапка: `{RESULTS_DIR}/{(row_index // 1000) * 1000:04d}/`
- Файл: `{apfid}.json` или `{apfid}__incorrect.json`

#### 8.4 `is_processed(row_index: int, apfid: str) -> bool`
- Если `ONLY_NOT_PROCESSED`: проверяет наличие `{apfid}.json` ИЛИ `{apfid}__incorrect.json`
- Иначе: `False`

#### 8.5 `process_experiment_group(args) -> list[dict]`

Функция-воркер для одного `experiment_id`. Принимает `(experiment_id, rows_list, global_indices)`.

Алгоритм:
1. Загрузить PDB: `load_structure(experiment_id)` → `structure`
2. Для каждой строки из `rows_list` (по порядку):
   a. Проверить `is_processed` → skip если нужно
   b. `chain = structure[0][chain_id]`
   c. `protein_residues = get_protein_residues(chain)`
   d. `chain_len = len(protein_residues)`
   e. `motif_residues = extract_motif_residues(chain, motif_start, motif_end)`
   f. `validate_motif(...)` → при ошибке → save `__incorrect.json` → continue
   g. `save_motif_to_pdb(motif_residues, chain_id, motif_pdb_path)` — только для STRIDE
   h. `save_chain_to_pdb(chain, chain_pdb_path)` — только для ProDy
   i. STRIDE: `run_stride(motif_pdb_path)` → `parse_stride_output` → `get_motif_ss_stats`
   j. SASA: `calc_sasa_isolated(motif_residues)` + `calc_sasa_in_chain(chain, motif_residue_ids)` → `surface_k`  
      *(передаём объекты из памяти, без чтения файлов)*
   k. Interactions: `calc_interactions(chain_pdb_path, motif_residue_ids)`
   l. `fasta = get_motif_fasta(motif_residues)`
   m. `rg = calc_rg(motif_residues)`
   n. `motif_weight = calc_motif_weight(motif_residues)`
   o. UniProt: `get_uniprot_data(uniprot_id)` → protein meta
   p. `protein_weight = calc_protein_weight(up_sequence)`
   q. `weight_percent = motif_weight / protein_weight * 100` (или `0.0` если `protein_weight == 0`)
   r. `gravy = calc_gravy(fasta)`
   s. Собрать итоговый `result_dict` со всеми полями
   t. Сохранить в JSON: `save_result(result_dict, row_index, apfid)`
   u. Удалить временные файлы: `motif_pdb_path`, `chain_pdb_path` (и `.stride` уже удалён в ЗАДАЧА 3.1)
3. При любом необработанном исключении:
   - Если `STOP_ON_ERROR`: перевозбудить
   - Иначе: залогировать в `errors.log` в формате:
     ```
     ==========================================================
     {apfid}
     {traceback}
     ```
   - Сохранить `__incorrect.json` с `apfid` + `error`

#### 8.6 Multiprocessing pool

```python
import multiprocessing, os

tasks = iter_experiment_groups(DATA_FILE, START_INDEX, STOP_INDEX)

with multiprocessing.Pool(processes=NUM_WORKERS or os.cpu_count()) as pool:
    for _ in pool.imap(process_experiment_group, tasks):
        pass  # результаты сохраняются в JSON внутри воркера; здесь просто дрейним итератор
```

- `pool.imap()` вместо `pool.map()` — генератор задач не материализуется целиком в памяти
- Каждый `experiment_id` = одна задача = один PDB-файл в одном процессе
- Порядок выполнения не гарантирован (это нормально)
- `NUM_WORKERS = None` → `os.cpu_count()`; каждый воркер держит в памяти не более одного PDB-файла

#### 8.7 Инициализация при старте
- Создать `RESULTS_DIR` и все подпапки если нет
- Создать `TEMP_DIR` если нет
- Вызвать `init_uniprot(UNIPROT_CACHE_DIR)`
- Настроить logging: **для каждого запуска создаётся отдельный лог-файл** с datetime в имени  
  - Формат datetime — Windows-совместимый (без символа `:`):  
    `processing_2026-05-12_14-30-05.log`  
  - Шаблон: `f"processing_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"`  
  - Дополнительно: вывод в stdout

---

### ЗАДАЧА 9: `aggregator.py` — JSON → CSV

**Файл:** `aggregator.py`  
**Зависит от:** ЗАДАЧА 1, ЗАДАЧА 8 (формат JSON)

> **Важно:** pandas не используем — работаем только через модуль `csv` для потоковой обработки.

#### 9.1 Логика агрегации

```python
# Параметры (из config):
START_INDEX = AGGREGATOR_START_INDEX
STOP_INDEX = AGGREGATOR_STOP_INDEX
```

1. Открываем входной CSV (`DATA_FILE`) через `csv.DictReader` и читаем построчно
2. Пропускаем строки до `START_INDEX`, останавливаемся на `STOP_INDEX`
3. Для каждой строки (по абсолютному индексу в файле):
   - Определяем путь JSON: `get_result_path(abs_index, apfid)`
   - **Если `{apfid}.json` существует**: читаем JSON, мержим поля поверх исходной строки (пересчитанные перезаписывают значения CSV), пишем в **correct**-файл
   - **Если `{apfid}__incorrect.json` существует**: читаем JSON, пишем `apfid` + `error` в **incorrect**-файл
   - **Если ни одного JSON нет**: пропускаем строку (не обработана)
4. Оба выходных файла открываются **заранее** (до цикла) через `csv.DictWriter` и пишутся потоково — строка за строкой, без накопления в памяти

#### 9.2 Формат имени выходных файлов
- Шаблон: `data/output_{start}_{stop}.csv` и `data/output_{start}_{stop}_incorrect.csv`
- `None` для `STOP_INDEX` → строка `"end"`: `output_0_end.csv`

#### 9.3 Заголовки выходных файлов
- **correct CSV**: все поля исходного CSV + поле `gravy` (если его нет в исходнике)
- **incorrect CSV**: только поля `apfid`, `error`
- Заголовки определяются по первой обработанной строке (как в примере `get_up_additionals.py`)

---

## Диаграмма зависимостей задач

```
ЗАДАЧА 0 (зависимости)
    │
    ├──► ЗАДАЧА 1 (config.py)
    │         │
    │    ┌────┴─────────────────────────────────────┐
    │    │                                           │
    ▼    ▼                                           │
ЗАДАЧА 2 (structure_tool.py)                        │
    │                                               │
    ├──► ЗАДАЧА 3 (stride_tool.py)                  │
    │         │                                     │
    ├──► ЗАДАЧА 4 (sasa_tool.py)                    │
    │         │                                     │
    ├──► ЗАДАЧА 5 (interactions_tool.py)            │
    │         │                                     │
    ├──► ЗАДАЧА 7 (gravy в structure_tool.py)       │
    │         │                                     │
    │    ЗАДАЧА 6 (uniprot_tool.py) ◄───────────────┘
    │         │
    └────┬────┘
         │
         ▼
   ЗАДАЧА 8 (processor.py)
         │
         ▼
   ЗАДАЧА 9 (aggregator.py)
```

---

## Технические уточнения

### PDB-путь
```python
def get_pdb_path(experiment_id: str) -> str:
    eid = experiment_id.lower()
    middle = eid[1:3]  # символы 2 и 3 (0-indexed: 1, 2)
    return f"{PDB_DIR}/{middle}/pdb{eid}.ent.gz"
```
Пример: `4E06` → `{PDB_DIR}/e0/pdb4e06.ent.gz`

### Временные файлы
- Мотив: `{TEMP_DIR}/{apfid}.pdb`
- Цепь: `{TEMP_DIR}/{apfid}_chain.pdb`
- Удаляются сразу после расчётов (в `finally`-блоке)

### Логирование ошибок
Для каждого запуска `processor.py` создаётся отдельный файл вида:
```
processing_2026-05-12_14-30-05.log
```
Формат datetime: `%Y-%m-%d_%H-%M-%S` — без двоеточий, Windows-совместимый.

Формат записи об ошибке внутри файла:
```
==========================================================
{apfid}
{полный traceback}
```

### Расчёт `motif_len`
`motif_len = motif_end - motif_start` (из CSV, не пересчитывается из реального кол-ва residue)

### Молекулярная масса с нестандартными остатками
- Стандартные (A–Y кроме B/Z) → BioPython `molecular_weight`
- `X` или `?` → 111.1 Да за каждый остаток (добавляется вручную к результату BioPython)
- Подход: заменить X/? на A в строке для BioPython + добавить `(111.1 - mass_of_A) * count_of_X`

### FreeSASA: версия API
Используем `freesasa.calcBioPDB(structure)` — принимает объект BioPython `Structure` напрямую из памяти, без записи в файл.  
Для суммирования SASA по нужным residue итерируемся по атомам результата `freesasa.Result`, фильтруя по `resnum`.

### ProDy: `Interactions`
```python
import prody
ag = prody.parsePDB(chain_pdb_path)
inter = prody.Interactions()
hb = inter.calcHydrogenBonds(ag)
hp = inter.calcHydrophobicContacts(ag)
ds = inter.calcDisulfideBonds(ag)
```
Каждый элемент результата содержит `resnum` обоих участников взаимодействия.

### GRAVY
```python
KD_SCALE = {'A': 1.8, 'R': -4.5, ...}  # Kyte-Doolittle
def calc_gravy(fasta: str) -> float:
    scores = [KD_SCALE[aa] for aa in fasta if aa in KD_SCALE]
    return sum(scores) / len(scores) if scores else 0.0
```

---

## Порядок реализации (рекомендуемый)

1. **ЗАДАЧА 0** — добавить зависимости
2. **ЗАДАЧА 1** — `config.py`
3. **ЗАДАЧА 2** — `structure_tool.py` (фундамент всего)
4. **ЗАДАЧА 7** — `calc_gravy` (простая, добавить в structure_tool.py)
5. **ЗАДАЧА 3** — `stride_tool.py`
6. **ЗАДАЧА 4** — `sasa_tool.py`
7. **ЗАДАЧА 5** — `interactions_tool.py`
8. **ЗАДАЧА 6** — `uniprot_tool.py`
9. **ЗАДАЧА 8** — `processor.py`
10. **ЗАДАЧА 9** — `aggregator.py`
