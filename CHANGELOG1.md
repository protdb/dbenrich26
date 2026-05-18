# CHANGELOG1 — Реализация проекта dbenrich26

**Дата:** 2026-05-12  
**Автор:** AI-агент по спецификации SPEC.md

---

## Общий прогресс

Выполнены **все 9 задач** из SPEC.md. Создано **8 Python-файлов** (1203 строки кода), обновлён `pyproject.toml`.

---

## ЗАДАЧА 0: Настройка зависимостей

**Файл:** `pyproject.toml` (обновлён)

Добавлены зависимости:
- `biopython (>=1.87,<2.0)` — загрузка/парсинг PDB, манипуляции с остатками, молекулярная масса
- `freesasa (>=2.2.1,<3.0.0)` — расчёт SASA напрямую из объектов BioPython
- `prody (>=2.6.1,<3.0.0)` — расчёт межостаточных взаимодействий (H-связи, гидрофобные контакты, дисульфидные мостики)

**Команда:** `poetry add biopython freesasa prody`  
**Залочены:** 5 пакетов (biopython 1.87, freesasa 2.2.1, prody 2.6.1, scipy 1.17.1, pyparsing 3.1.1)

---

## ЗАДАЧА 1: `config.py`

**Файл:** `config.py` — 25 строк

Вынесены все константы:
- **Пути:** `PDB_DIR`, `DATA_FILE`, `RESULTS_DIR`, `TEMP_DIR`, `UNIPROT_CACHE_DIR`
- **Управление:** `START_INDEX`, `STOP_INDEX`, `ONLY_NOT_PROCESSED`, `STOP_ON_ERROR`
- **Параметры:** `NUM_WORKERS`, `MISSING_RESIDUE_THRESHOLD` (10%), `COIL_TAIL_MAX_LEN` (3)
- **Агрегатор:** `AGGREGATOR_START_INDEX`, `AGGREGATOR_STOP_INDEX`

---

## ЗАДАЧА 2: `structure_tool.py` + ЗАДАЧА 7: GRAVY

**Файл:** `structure_tool.py` — 235 строк

### Функции:

| Функция | Описание |
|---------|----------|
| `get_pdb_path(experiment_id)` | Конструирует путь в архиве PDB (шаблон: `{PDB_DIR}/{mid}/pdb{eid}.ent.gz`) |
| `load_structure(experiment_id)` | Загружает структуру через `Bio.PDB.PDBParser` из `.ent.gz` |
| `is_protein_residue(residue)` | Фильтр: стандартные + модифицированные (с CA), исключая воду/лиганды |
| `get_protein_residues(chain)` | Все белковые остатки цепи |
| `get_chain_len(chain)` | Длина цепи (только белковые остатки) |
| `extract_motif_residues(chain, start, end)` | Вырезает остатки мотива по sequence numbers |
| `validate_motif(residues, start, end)` | Проверка: не пустой, пропуски ≤ 10% |
| `save_motif_to_pdb(residues, chain_id, path)` | Сохраняет мотив в PDB-файл (для STRIDE) |
| `save_chain_to_pdb(chain, path)` | Сохраняет цепь в PDB-файл (для ProDy) |
| `calc_rg(residues)` | Радиус гирации по Cα-атомам |
| `get_motif_fasta(residues)` | Однобуквенная последовательность, нестандартные → `?` |
| `calc_motif_weight(residues)` | Молекулярная масса мотива (BioPython + коррекция на `?` = 111.1 Да) |
| `calc_protein_weight(sequence)` | Масса из строки FASTA |
| `calc_gravy(fasta)` | **GRAVY** по шкале Kyte-Doolittle (`?` / `X` пропускаются) |

### Классы исключений:
- `MotifValidationError` — пустой мотив или превышен лимит пропусков

---

## ЗАДАЧА 3: `stride_tool.py`

**Файл:** `stride_tool.py` — 150 строк

### Функции:

| Функция | Описание |
|---------|----------|
| `run_stride(pdb_path)` | Запускает бинарник `./stride`, возвращает raw-текст. `.stride`-файл удаляется |
| `parse_stride_output(raw)` | Парсит ASG-строки, ремаппит T/B/G → C (coil) |
| `get_secondary_structure_elements(parsed)` | Схлопывает остатки в непрерывные элементы SS |
| `get_motif_ss_stats(elements, start, end, max_coil_tail)` | Фильтрует элементы по границам мотива, отсекает хвостовые coil'ы (≤3), считает статистику |

**Логика хвостов:** первый и последний элемент типа C длиной ≤3 не учитываются (однократно, не рекурсивно).

---

## ЗАДАЧА 4: `sasa_tool.py`

**Файл:** `sasa_tool.py` — 66 строк

### Функции:

| Функция | Описание |
|---------|----------|
| `calc_sasa_isolated(residues)` | SASA мотива в изоляции через `freesasa.calcBioPDB()` — строго из памяти, без записи на диск |
| `calc_sasa_in_chain(chain, residue_ids)` | SASA полной цепи, суммирование только по атомам остатков мотива |
| `calc_surface_k(isolated, in_chain)` | `surface_k = in_chain / isolated`, при `isolated == 0` → 0.0 |

**Внутренняя утилита:** `_build_single_chain_structure(residues)` — собирает временную in-memory структуру.

---

## ЗАДАЧА 5: `interactions_tool.py`

**Файл:** `interactions_tool.py` — 74 строки

### Функции:

| Функция | Описание |
|---------|----------|
| `calc_interactions(chain_pdb_path, motif_ids)` | Загружает цепь через `prody.parsePDB`, вычисляет три типа взаимодействий и разбивает на inner/outer |

**Типы:** H-связи (`calcHydrogenBonds`), гидрофобные (`calcHydrophobicContacts`), дисульфидные (`calcDisulfideBonds`).  
**Критерии:** inner = оба resnum в мотиве; outer = один в мотиве; оба снаружи — игнорируются.

---

## ЗАДАЧА 6: `uniprot_tool.py`

**Файл:** `uniprot_tool.py` — 156 строк → 90 строк (после рефакторинга)

### Функции (оригинал):

| Функция | Описание |
|---------|----------|
| `init_uniprot(cache_dir)` | Инициализирует кэш UniProt |
| `get_uniprot_data(uniprot_id)` | Извлекает: protein_name, GO-термы, organism, sequence |
| `enrich_with_partners_and_sites(data, exp_id, chain_id, uniprot_id, start, end)` | Дополняет: partner proteins (через PDB API), binding sites (через `get_position_features`) |

**Рецепт партнёров:** перенесён из `examples/get_up_additionals.py` (строки 68–77) — `jmespath`-запрос `DATA_PLAIN` к PDB API, дедупликация по UniProt ID.

**Fallback:** при недоступном UniProt — исходные поля сохраняются, protein_weight = 0.

---

## ЗАДАЧА 8: `processor.py`

**Файл:** `processor.py` — 375 строк (самый большой файл)

### Ключевые компоненты:

1. **`iter_experiment_groups(filepath, start, stop)`** — потоковый генератор на `csv.DictReader`.
   - Читает строки, группирует по `experiment_id`
   - Применяет `START_INDEX`/`STOP_INDEX` без загрузки всего файла
   - Выдаёт `(experiment_id, [rows], [indices])` — одну группу за раз

2. **`get_result_path()`** — путь к JSON в подпапках `0000/`, `1000/`, ...
3. **`is_processed()`** — проверка наличия JSON (если `ONLY_NOT_PROCESSED`)
4. **`process_experiment_group(args)`** — воркер для одного PDB-файла:
   - Загружает PDB один раз (все мотивы группы)
   - Для каждого мотива: save PDB → STRIDE → SASA → ProDy → Rg/масса → GRAVY → UniProt
   - Сохраняет JSON или `__incorrect.json`
   - Удаляет временные файлы в `finally`
5. **`main()`** — инициализация директорий, кэша, лога, пул `multiprocessing.Pool.imap()`

### Логирование:
- **Формат файла:** `processing_2026-05-12_22-50-00.log` (datetime без двоеточий)
- **Формат ошибки:** `===== {apfid} \n {traceback}`
- **Два хендлера:** файл + stdout

---

## ЗАДАЧА 9: `aggregator.py`

**Файл:** `aggregator.py` — 122 строки

### Алгоритм:

1. Читает `data_SORTED.csv` через `csv.DictReader` построчно
2. Для каждой строки (в срезе `AGGREGATOR_START..STOP`):
   - Находит JSON по `get_result_path(abs_index, apfid)`
   - `{apfid}.json` → мержит поля поверх исходной строки → пишет в **correct** CSV
   - `{apfid}__incorrect.json` → пишет `apfid` + `error` в **incorrect** CSV
   - Нет JSON → пропуск (не обработана)
3. Потоковая запись: `csv.DictWriter`, без накопления в памяти
4. **Выходные файлы:** `data/output_{start}_{stop}.csv` и `data/output_{start}_{stop}_incorrect.csv`

---

## Итого по проекту

```
dbenrich26/
├── pyproject.toml          # 23 стр. — зависимости (biopython, freesasa, prody добавлены)
├── config.py               # 25 стр. — конфигурация
├── structure_tool.py       # 235 стр. — ядро: PDB, residue, Rg, масса, FASTA, GRAVY
├── stride_tool.py          # 150 стр. — STRIDE: запуск, парсинг, статистика
├── sasa_tool.py            # 66 стр. — FreeSASA из памяти
├── interactions_tool.py    # 74 стр. — ProDy: inner/outer взаимодействия
├── uniprot_tool.py         # 156 стр. — UniProt: метаданные, партнёры, binding sites
├── processor.py            # 375 стр. — основной multiprocessing pipeline
├── aggregator.py           # 122 стр. — JSON → CSV агрегатор
├── stride                  # бинарник STRIDE (без изменений)
└── data/
    ├── data_SORTED.csv     # 1.2M строк (без изменений)
    └── results/            # создаётся автоматически
```

**Всего:** 1137 строк Python-кода, 8 новых файлов.

---

## Исправление: удаление PDB API из `uniprot_tool.py`

**Дата:** 2026-05-18  
**Причина:** `ModuleNotFoundError: No module named 'get_pdb_api_data'` — модуль `get_pdb_api_data` (и `apfid`) не существовал в проекте.

### Что изменилось

#### `uniprot_tool.py` — 156 → 90 строк (-66 строк)

**Удалено:**
- `from get_pdb_api_data import get_single` — несуществующий модуль
- `from apfid import Apfid` — несуществующий модуль
- `import jmespath` — нужен был только для PDB API
- Константа `DATA_PLAIN` — JMESPath-запрос к PDB API
- Функция `enrich_with_partners_and_sites()` — 170 строк, вся логика получения партнёров через `Apfid` + `get_single` + `jmespath`

**Добавлено:**
- Функция `enrich_binding_sites(up_data, uniprot_id, motif_start, motif_end)`:
  - Только сайты связывания через `UniprotData.get_position_features()`
  - `protein_partners` всегда устанавливается в `''` (PDB API не используется)
  - Сигнатура проще: не требует `experiment_id` и `chain_id`

#### `processor.py` — косметические правки
- Импорт: `enrich_with_partners_and_sites` → `enrich_binding_sites`
- Вызов: убраны параметры `experiment_id, chain_id`

**Не затронуто:** логика сайтов связывания (`has_binding_sites`, `binding_sites`) — осталась без изменений, используется тот же `UniprotData.get_position_features()`.

### Изменение метрик
| Параметр | До | После |
|---|---|---|
| `uniprot_tool.py` | 156 строк | 90 строк |
| Всего кода | 1203 строки | 1137 строк |
| Файлов | 8 | 8 |
