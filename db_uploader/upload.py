#!/usr/bin/env python3
"""Upload CSV output from aggregator.py into public.motif table in PostgreSQL.

Usage (Docker):
    docker run --rm -e PG_DSN="postgresql://...:5432/db" -v /path/to/csv:/data glucksistemi/psskb_uploader

Environment:
    PG_DSN      — PostgreSQL connection string (required)
    BATCH_SIZE  — rows between commits (default: 500)
"""

import csv
import glob
import math
import os
import shutil
import sys

import psycopg2
from tabulate import tabulate

PROGRESS_EVERY = 100
DONE_DIR = '/data/done'
BATCH_SIZE = 500  # default, overridden from env in main()
WRITE_MODE = 'update'  # default, overridden from env in main()

REAL_FIELDS = {
    'resolution', 'rg', 'sasa', 'surface_k',
    'motif_weight', 'protein_weight', 'weight_percent', 'gravy',
}

INT_FIELDS = {
    'model_id', 'chain_len', 'motif_len', 'motif_start', 'motif_end',
    'inner_hb', 'inner_hp', 'inner_ds',
    'outer_hb', 'outer_hp', 'outer_ds',
    'helix_count', 'strand_count', 'coil_count',
    'helix_length', 'strand_length', 'coil_length',
}

BOOL_FIELDS = {
    'full_chain', 'single_chain',
}

# Threshold checks: field → max value (exclusive, skipped if value > max)
# NaN/None values always pass the check.
FIELD_LIMITS = {
    'inner_hb': 78,
    'inner_hp': 93,
    'helix_count': 14,
    'strand_count': 16,
}


def check_field_limits(row: dict) -> bool:
    """Return True if all FIELD_LIMITS checks pass for the given row."""
    for field, limit in FIELD_LIMITS.items():
        val = row.get(field, '')
        if val == '' or val is None or (isinstance(val, str) and val.upper() == 'NULL'):
            continue
        try:
            if float(val) > limit:
                return False
        except (ValueError, TypeError):
            continue
    return True


def convert_field(field_name: str, value: str):
    """Convert CSV string to Python type matching the DB column type.

    - real fields: '' or 'NULL' → NaN (not NULL — PG real has no NULL)
    - int fields: '' or 'NULL' → None (SQL NULL)
    - bool fields: '' or 'NULL' → None
    - text/varchar fields: 'NULL' → None, '' stays ''
    """
    is_null = value is None or (isinstance(value, str) and value.upper() == 'NULL')

    if field_name in REAL_FIELDS:
        if is_null or value == '':
            return float('nan')
        return float(value)

    if field_name in INT_FIELDS:
        if is_null or value == '':
            return None
        return int(value)

    if field_name in BOOL_FIELDS:
        if is_null or value == '':
            return None
        return value.lower() in ('true', 't', '1', 'yes')

    # Text/varchar: NULL string → None, empty → '' (kept for NOT NULL columns)
    if is_null:
        return None
    return value


def count_data_lines(filepath: str) -> int:
    """Count data rows in a CSV file (total lines minus header), streaming."""
    count = 0
    with open(filepath, 'r', encoding='utf-8') as f:
        for _ in f:
            count += 1
    return max(count - 1, 0)


def detect_file_type(filename: str) -> str:
    """Return 'incorrect' if filename contains '_incorrect', else 'correct'."""
    return 'incorrect' if '_incorrect' in filename else 'correct'


def build_update_sql(fieldnames: list[str]):
    """Build UPDATE SQL once per file.

    Returns (sql_string, set_fields_list).
    Excludes 'id' and 'apfid' from the SET clause.
    """
    set_fields = [f for f in fieldnames if f.lower() not in ('id', 'apfid')]
    set_clause = ', '.join(f'{f} = %s' for f in set_fields)
    sql = f'UPDATE public.motif SET {set_clause} WHERE apfid = %s'
    return sql, set_fields


def make_update_params(row: dict, field_list: list[str]):
    """Build a parameter tuple for the UPDATE statement.

    Uses convert_field() to cast each value to the correct Python type
    matching the DB column (real → NaN for empty, int → None, etc.).
    """
    values = [convert_field(f, row.get(f, '')) for f in field_list]
    values.append(row.get('apfid', ''))
    return tuple(values)


def build_insert_sql(fieldnames: list[str]) -> tuple[str, list[str]]:
    """Build INSERT SQL once per file.

    Returns (sql_string, insert_fields_list).
    Excludes 'id' from the column list (autoincrement).
    'apfid' IS included.
    """
    insert_fields = [f for f in fieldnames if f.lower() != 'id']
    columns = ', '.join(insert_fields)
    placeholders = ', '.join('%s' for _ in insert_fields)
    sql = f'INSERT INTO public.motif ({columns}) VALUES ({placeholders})'
    return sql, insert_fields


def make_insert_params(row: dict, field_list: list[str]) -> tuple:
    """Build parameter tuple for INSERT. Uses convert_field() for type casting."""
    return tuple(convert_field(f, row.get(f, '')) for f in field_list)


def fetch_existing_motif_types(conn, apfids: list[str]) -> dict[str, str]:
    """Batch SELECT: returns {apfid: motif_type} for all found records.
    Records not in DB are NOT included in the dict.
    """
    if not apfids:
        return {}
    with conn.cursor() as cur:
        cur.execute(
            'SELECT apfid, motif_type FROM public.motif WHERE apfid = ANY(%s)',
            (apfids,)
        )
        return {row[0]: row[1] for row in cur.fetchall()}


def collect_apfids(filepath: str) -> list[str]:
    """First pass: collect all apfid values from a correct CSV file."""
    apfids = []
    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            apfid = row.get('apfid', '').strip()
            if apfid:
                apfids.append(apfid)
    return apfids


def process_file(conn, write_mode: str, filepath: str) -> dict:
    """Process a single CSV file: UPDATE for correct, DELETE for incorrect.

    Returns a stats dict:
        {filename, type, total, ok, not_found, errors}
    """
    filename = os.path.basename(filepath)
    total_lines = count_data_lines(filepath)
    file_type = detect_file_type(filename)

    stats = {
        'filename': filename,
        'type': file_type,
        'total': total_lines,
        'ok': 0,
        'not_found': 0,
        'skipped': 0,
        'skipped_protected': 0,
        'errors': 0,
    }

    if total_lines == 0:
        print(f"[{filename}] пустой файл, пропускаем")
        shutil.move(filepath, os.path.join(DONE_DIR, filename))
        print(f"[{filename}] → перемещён в done/{filename}")
        return stats

    # In upsert mode: skip incorrect files entirely
    if write_mode == 'upsert' and file_type == 'incorrect':
        print(f"[{filename}] incorrect-файл, режим upsert — пропускаем (без DELETE)")
        shutil.move(filepath, os.path.join(DONE_DIR, filename))
        print(f"[{filename}] → перемещён в done/{filename}")
        return stats

    # Upsert mode: first pass to collect all apfids, then batch query DB
    if write_mode == 'upsert' and file_type == 'correct':
        apfids = collect_apfids(filepath)
        existing_motif_types = fetch_existing_motif_types(conn, apfids)
        print(f"[{filename}] найдено в БД: {len(existing_motif_types)} / {len(apfids)} записей")

    # SQL templates
    delete_sql = 'DELETE FROM public.motif WHERE apfid = %s'
    update_sql = None
    update_fields = None

    # Skip file – lazily opened on first skipped row
    skip_path = os.path.join(DONE_DIR, f'{os.path.splitext(filename)[0]}_skip.csv')
    skip_fh = None
    skip_writer = None

    skip_protected_path = os.path.join(DONE_DIR, f'{os.path.splitext(filename)[0]}_protected.csv')
    skip_protected_fh = None
    skip_protected_writer = None

    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        csv_fieldnames = reader.fieldnames

        if file_type == 'correct':
            update_sql, update_fields = build_update_sql(reader.fieldnames)
            if write_mode == 'upsert':
                insert_sql, insert_fields = build_insert_sql(reader.fieldnames)

        for line_no, row in enumerate(reader, start=1):
            # Progress every PROGRESS_EVERY rows and on the last row
            if line_no % PROGRESS_EVERY == 0 or line_no == total_lines:
                print(f"[{filename}] строка {line_no} / {total_lines}")

            # Field limit checks – skip row if exceeded
            if not check_field_limits(row):
                stats['skipped'] += 1
                if skip_writer is None:
                    skip_fh = open(skip_path, 'w', newline='', encoding='utf-8')
                    skip_writer = csv.DictWriter(skip_fh, fieldnames=csv_fieldnames)
                    skip_writer.writeheader()
                skip_writer.writerow(row)
                continue

            try:
                with conn.cursor() as cur:
                    if write_mode == 'upsert':
                        apfid = row.get('apfid', '').strip()
                        if apfid not in existing_motif_types:
                            # INSERT
                            params = make_insert_params(row, insert_fields)
                            cur.execute(insert_sql, params)
                            stats['ok'] += 1
                        elif existing_motif_types[apfid] == 'Full chain':
                            # UPDATE
                            params = make_update_params(row, update_fields)
                            cur.execute(update_sql, params)
                            if cur.rowcount == 0:
                                print(f"  WARNING: [{filename}] строка {line_no}: apfid не найден: {apfid}")
                                stats['not_found'] += 1
                            else:
                                stats['ok'] += 1
                        else:
                            # Protected: motif_type != 'Full chain'
                            stats['skipped_protected'] += 1
                            if skip_protected_writer is None:
                                skip_protected_fh = open(skip_protected_path, 'w', newline='', encoding='utf-8')
                                skip_protected_writer = csv.DictWriter(skip_protected_fh, fieldnames=csv_fieldnames)
                                skip_protected_writer.writeheader()
                            skip_protected_writer.writerow(row)
                    else:
                        # Original update mode
                        if file_type == 'correct':
                            params = make_update_params(row, update_fields)
                            cur.execute(update_sql, params)
                        else:
                            cur.execute(delete_sql, (row.get('apfid', ''),))

                        if cur.rowcount == 0:
                            apfid = row.get('apfid', '?')
                            print(f"  WARNING: [{filename}] строка {line_no}: apfid не найден: {apfid}")
                            stats['not_found'] += 1
                        else:
                            stats['ok'] += 1

            except Exception as e:
                print(f"  ERROR: [{filename}] строка {line_no}: {e}")
                conn.rollback()
                stats['errors'] += 1

            # Periodic commit
            if line_no % BATCH_SIZE == 0:
                conn.commit()

    # Final commit for the file
    conn.commit()

    # Close skip-protected file if we wrote any
    if skip_protected_fh is not None:
        skip_protected_fh.close()
        print(f"[{filename}] защищено {stats['skipped_protected']} строк → {os.path.basename(skip_protected_path)}")

    # Close skip file if we wrote any
    if skip_fh is not None:
        skip_fh.close()
        print(f"[{filename}] пропущено {stats['skipped']} строк → {os.path.basename(skip_path)}")

    # Move file to done/ only after successful processing (including final commit)
    shutil.move(filepath, os.path.join(DONE_DIR, filename))
    print(f"[{filename}] готово: обработано {stats['ok'] + stats['not_found']}: "
          f"OK {stats['ok']}, не найдено {stats['not_found']}, "
          f"пропущено {stats['skipped']}, защищено {stats['skipped_protected']}, "
          f"ошибок {stats['errors']}")
    print(f"[{filename}] → перемещён в done/{filename}")

    return stats


def main():
    global BATCH_SIZE, WRITE_MODE

    pg_dsn = os.environ.get('PG_DSN')
    if not pg_dsn:
        print("ERROR: переменная PG_DSN не задана", file=sys.stderr)
        sys.exit(1)

    try:
        BATCH_SIZE = int(os.environ.get('BATCH_SIZE', str(BATCH_SIZE)))
    except ValueError:
        print(f"ERROR: BATCH_SIZE должен быть целым числом, получено: {os.environ.get('BATCH_SIZE')}",
              file=sys.stderr)
        sys.exit(1)

    raw_mode = os.environ.get('WRITE_MODE', 'update').lower().strip()
    if raw_mode not in ('update', 'upsert'):
        print(f"ERROR: WRITE_MODE должен быть 'update' или 'upsert', получено: {os.environ.get('WRITE_MODE')}",
              file=sys.stderr)
        sys.exit(1)
    WRITE_MODE = raw_mode

    files = sorted(glob.glob('/data/output*.csv'))
    files = [f for f in files if os.path.isfile(f)]

    if not files:
        print("Нет CSV-файлов в /data")
        sys.exit(0)

    os.makedirs(DONE_DIR, exist_ok=True)

    print(f"Подключение к БД...")
    conn = psycopg2.connect(pg_dsn)
    print(f"OK. batch_size={BATCH_SIZE}, write_mode={WRITE_MODE}, файлов: {len(files)}")
    print()

    all_stats = []

    try:
        for filepath in files:
            filename = os.path.basename(filepath)
            print(f"--- {filename} ---")
            stats = process_file(conn, WRITE_MODE, filepath)
            all_stats.append(stats)
            print()
    finally:
        conn.close()

    # --- Final summary table ---
    print("=" * 80)
    print("СВОДНАЯ СТАТИСТИКА")
    print("=" * 80)

    table_rows = []
    sum_total = sum_ok = sum_not_found = sum_skipped = sum_skipped_protected = sum_errors = 0

    for s in all_stats:
        status = '✓' if s['errors'] == 0 else '✗'
        table_rows.append([
            s['filename'],
            s['type'],
            s['total'],
            s['ok'],
            s['not_found'],
            s['skipped'],
            s['skipped_protected'],
            s['errors'],
            status,
        ])
        sum_total += s['total']
        sum_ok += s['ok']
        sum_not_found += s['not_found']
        sum_skipped += s['skipped']
        sum_skipped_protected += s['skipped_protected']
        sum_errors += s['errors']

    all_ok = sum_errors == 0
    table_rows.append([
        '**SUM**', '',
        sum_total, sum_ok, sum_not_found, sum_skipped, sum_skipped_protected, sum_errors,
        '✓' if all_ok else '✗',
    ])

    headers = ['Файл', 'Тип', 'Строк', 'OK', 'Не найдено', 'Пропущено', 'Защищено', 'Ошибок', 'Статус']
    print(tabulate(table_rows, headers=headers, tablefmt='pipe'))
    print()

    sys.exit(0 if all_ok else 1)


if __name__ == '__main__':
    main()
