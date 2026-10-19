"""Aggregator: converts JSON results from processor.py into final CSV files."""
import csv
import json
import os
import sys
from typing import Optional

import config as cfg


def get_result_path(row_index: int, apfid: str, incorrect: bool = False) -> str:
    """Возвращает путь к JSON-файлу результата (аналогично processor.py)."""
    subdir_num = (row_index // 1000) * 1000
    subdir = os.path.join(cfg.RESULTS_DIR, f"{subdir_num:04d}")
    suffix = "__incorrect" if incorrect else ""
    return os.path.join(subdir, f"{apfid}{suffix}.json")


def _format_stop(stop: Optional[int]) -> str:
    return str(stop) if stop is not None else "end"


def _get_fieldnames(row: dict, result: dict) -> list:
    """
    Определяет заголовки для выходного CSV: все поля исходной строки + поле gravy,
    если его нет в исходнике.
    """
    fieldnames = list(row.keys())
    if 'gravy' not in fieldnames:
        fieldnames.append('gravy')
    return fieldnames


def main():
    start = cfg.AGGREGATOR_START_INDEX
    stop = cfg.AGGREGATOR_STOP_INDEX
    stop_str = _format_stop(stop)

    output_correct = os.path.join("data", f"output_{start}_{stop_str}.csv")
    output_incorrect = os.path.join("data", f"output_{start}_{stop_str}_incorrect.csv")

    correct_writer = None
    incorrect_writer = None
    correct_file = None
    incorrect_file = None

    correct_count = 0
    incorrect_count = 0
    skipped_count = 0

    print(f"Aggregating: START={start}, STOP={stop_str}")
    print(f"Correct   → {output_correct}")
    print(f"Incorrect → {output_incorrect}")

    with open(cfg.DATA_FILE, newline='', encoding='utf-8') as infile:
        reader = csv.DictReader(infile)

        for abs_index, row in enumerate(reader):
            if abs_index < start:
                continue
            if stop is not None and abs_index >= stop:
                break

            apfid = row['apfid']

            # Пробуем correct JSON
            correct_path = get_result_path(abs_index, apfid, incorrect=False)
            incorrect_path = get_result_path(abs_index, apfid, incorrect=True)

            if os.path.exists(correct_path):
                # Читаем JSON и мержим поверх исходной строки
                with open(correct_path, 'r', encoding='utf-8') as jf:
                    result = json.load(jf)

                # Сливаем: пересчитанные поля перезаписывают исходные
                merged = dict(row)
                merged.update(result)

                # Открываем correct CSV при первой записи
                if correct_writer is None:
                    fieldnames = _get_fieldnames(row, result)
                    correct_file = open(output_correct, 'w', newline='', encoding='utf-8')
                    correct_writer = csv.DictWriter(correct_file, fieldnames=fieldnames)
                    correct_writer.writeheader()

                correct_writer.writerow(merged)
                correct_count += 1

            elif os.path.exists(incorrect_path):
                # Читаем incorrect JSON и пишем apfid + error
                with open(incorrect_path, 'r', encoding='utf-8') as jf:
                    err_data = json.load(jf)

                if incorrect_writer is None:
                    incorrect_file = open(output_incorrect, 'w', newline='', encoding='utf-8')
                    incorrect_writer = csv.DictWriter(
                        incorrect_file, fieldnames=['apfid', 'error']
                    )
                    incorrect_writer.writeheader()

                incorrect_writer.writerow({
                    'apfid': err_data.get('apfid', apfid),
                    'error': err_data.get('error', ''),
                })
                incorrect_count += 1

            else:
                print(correct_path)
                print(incorrect_path)
                skipped_count += 1

    # Закрываем файлы
    if correct_file:
        correct_file.close()
    if incorrect_file:
        incorrect_file.close()

    print(f"\nDone.")
    print(f"  Correct rows:   {correct_count}")
    print(f"  Incorrect rows: {incorrect_count}")
    print(f"  Skipped (no result): {skipped_count}")


if __name__ == "__main__":
    main()