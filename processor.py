"""Main multiprocessing pipeline for processing PDB motifs."""
import csv
import json
import logging
import multiprocessing
import os
import sys
import traceback
from datetime import datetime
from typing import Iterator, Tuple, List, Optional, Dict

import config as cfg
from structure_tool import (
    load_structure, get_protein_residues, get_chain_len,
    extract_motif_residues, validate_motif, MotifValidationError,
    save_motif_to_pdb, save_chain_to_pdb,
    calc_rg, calc_motif_weight, calc_protein_weight,
    get_motif_fasta, calc_gravy,
)
from stride_tool import run_stride, parse_stride_output, get_secondary_structure_elements, get_motif_ss_stats
from sasa_tool import calc_sasa_isolated, calc_sasa_in_chain, calc_surface_k
from interactions_tool import calc_interactions
from uniprot_tool import init_uniprot, get_uniprot_data, enrich_with_partners_and_sites

# ─── Логирование ─────────────────────────────────────────────────────────────


def setup_logging() -> str:
    """Настраивает логирование: файл с datetime + stdout. Возвращает имя лог-файла."""
    log_name = f"processing_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_name, encoding='utf-8'),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return log_name


def log_error(apfid: str, exc_info) -> None:
    """Логирует ошибку в формате: === apfid === traceback."""
    logging.error(
        "==========================================================\n"
        "%s\n%s",
        apfid,
        "".join(traceback.format_exception(*exc_info)),
    )


# ─── Чтение входных данных (генератор групп) ──────────────────────────────


def iter_experiment_groups(
    filepath: str, start_index: int = 0, stop_index: Optional[int] = None,
) -> Iterator[Tuple[str, List[Dict], List[int]]]:
    """
    Читает CSV построчно, группируя строки по experiment_id.
    Файл должен быть отсортирован по experiment_id.
    Выдаёт (experiment_id, rows_list, global_indices_list) — по одной группе за раз.
    """
    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        current_group_id = None
        current_rows: List[Dict] = []
        current_indices: List[int] = []

        for abs_index, row in enumerate(reader):
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


# ─── Пути результатов ──────────────────────────────────────────────────────


def get_result_path(row_index: int, apfid: str, incorrect: bool = False) -> str:
    """Возвращает путь к JSON-файлу результата для данного индекса."""
    subdir_num = (row_index // 1000) * 1000
    subdir = os.path.join(cfg.RESULTS_DIR, f"{subdir_num:04d}")
    os.makedirs(subdir, exist_ok=True)
    suffix = "__incorrect" if incorrect else ""
    return os.path.join(subdir, f"{apfid}{suffix}.json")


def is_processed(row_index: int, apfid: str) -> bool:
    """Проверяет, существует ли уже результат для данного apfid."""
    if not cfg.ONLY_NOT_PROCESSED:
        return False
    correct_path = get_result_path(row_index, apfid, incorrect=False)
    incorrect_path = get_result_path(row_index, apfid, incorrect=True)
    return os.path.exists(correct_path) or os.path.exists(incorrect_path)


def save_result(result_dict: dict, row_index: int, apfid: str) -> None:
    """Сохраняет словарь результата в JSON-файл."""
    path = get_result_path(row_index, apfid, incorrect=False)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(result_dict, f, ensure_ascii=False)


def save_incorrect_result(apfid: str, error_msg: str, row_index: int) -> None:
    """Сохраняет JSON для некорректной структуры."""
    path = get_result_path(row_index, apfid, incorrect=True)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump({"apfid": apfid, "error": error_msg}, f, ensure_ascii=False)


# ─── Воркер ────────────────────────────────────────────────────────────────


def process_experiment_group(args: Tuple[str, List[Dict], List[int]]) -> None:
    """
    Обрабатывает один experiment_id: загружает PDB, вырезает все мотивы,
    выполняет все расчёты, сохраняет JSON.
    Вызывается внутри multiprocessing.Pool.
    """
    experiment_id, rows, global_indices = args

    # Загружаем PDB один раз для всех мотивов этого experiment_id
    try:
        structure = load_structure(experiment_id)
    except FileNotFoundError as e:
        for row, idx in zip(rows, global_indices):
            apfid = row['apfid']
            err = f"PDB file not found: {e}"
            save_incorrect_result(apfid, err, idx)
            if cfg.STOP_ON_ERROR:
                raise
            log_error(apfid, sys.exc_info())
        return
    except Exception as e:
        for row, idx in zip(rows, global_indices):
            apfid = row['apfid']
            err = f"Failed to load PDB: {e}"
            save_incorrect_result(apfid, err, idx)
            if cfg.STOP_ON_ERROR:
                raise
            log_error(apfid, sys.exc_info())
        return

    model = structure[0]

    for row, global_idx in zip(rows, global_indices):
        apfid = row['apfid']
        chain_id = row['chain_id']

        # Проверка на уже обработанное
        if is_processed(global_idx, apfid):
            logging.info(f"SKIP (already processed): {apfid}")
            continue

        motif_start = int(row['motif_start'])
        motif_end = int(row['motif_end'])
        uniprot_id = row.get('uniprot_id', '')

        # Получаем цепь
        try:
            chain = model[chain_id]
        except KeyError:
            err = f"Chain {chain_id} not found in structure {experiment_id}"
            save_incorrect_result(apfid, err, global_idx)
            log_error(apfid, sys.exc_info())
            continue

        # Белковые остатки цепи
        try:
            protein_residues = get_protein_residues(chain)
            chain_len = len(protein_residues)
        except Exception as e:
            err = f"Failed to get protein residues: {e}"
            save_incorrect_result(apfid, err, global_idx)
            log_error(apfid, sys.exc_info())
            continue

        # Вырезаем мотив
        try:
            motif_residues = extract_motif_residues(chain, motif_start, motif_end)
            validate_motif(motif_residues, motif_start, motif_end)
        except MotifValidationError as e:
            err = f"Motif validation failed: {e}"
            save_incorrect_result(apfid, err, global_idx)
            log_error(apfid, sys.exc_info())
            continue
        except Exception as e:
            err = f"Failed to extract motif: {e}"
            save_incorrect_result(apfid, err, global_idx)
            log_error(apfid, sys.exc_info())
            continue

        # Пути для временных файлов
        motif_pdb_path = os.path.join(cfg.TEMP_DIR, f"{apfid}.pdb")
        chain_pdb_path = os.path.join(cfg.TEMP_DIR, f"{apfid}_chain.pdb")
        os.makedirs(cfg.TEMP_DIR, exist_ok=True)

        try:
            # Сохраняем PDB-фрагменты (нужны для STRIDE и ProDy)
            save_motif_to_pdb(motif_residues, chain_id, motif_pdb_path)
            save_chain_to_pdb(chain, chain_pdb_path)

            # 1. STRIDE — вторичная структура
            stride_raw = run_stride(motif_pdb_path)
            ss_stats = {}
            if stride_raw:
                parsed = parse_stride_output(stride_raw)
                if parsed:
                    elements = get_secondary_structure_elements(parsed)
                    ss_stats = get_motif_ss_stats(elements, motif_start, motif_end)
            else:
                ss_stats = {
                    "helix_count": 0, "strand_count": 0, "coil_count": 0,
                    "helix_length": 0, "strand_length": 0, "coil_length": 0,
                }

            # 2. SASA
            try:
                sasa_isolated = calc_sasa_isolated(motif_residues)
                motif_residue_ids = [r.id[1] for r in motif_residues]
                sasa_in_chain = calc_sasa_in_chain(chain, motif_residue_ids)
                surface_k = calc_surface_k(sasa_isolated, sasa_in_chain)
            except Exception as e:
                logging.warning(f"SASA calculation failed for {apfid}: {e}")
                sasa_isolated = 0.0
                surface_k = 0.0

            # 3. Interactions — ProDy
            try:
                motif_residue_ids_set = set(r.id[1] for r in motif_residues)
                interactions = calc_interactions(chain_pdb_path, motif_residue_ids_set)
            except Exception as e:
                logging.warning(f"Interactions calculation failed for {apfid}: {e}")
                interactions = {
                    "inner_hb": 0, "inner_hp": 0, "inner_ds": 0,
                    "outer_hb": 0, "outer_hp": 0, "outer_ds": 0,
                }

            # 4. FASTA, Rg, масса
            fasta = get_motif_fasta(motif_residues)
            rg = calc_rg(motif_residues)
            motif_weight = calc_motif_weight(motif_residues)

            # 5. GRAVY
            gravy = calc_gravy(fasta)

            # 6. UniProt
            up_data = get_uniprot_data(uniprot_id)
            if up_data is not None:
                try:
                    up_data = enrich_with_partners_and_sites(
                        up_data, experiment_id, chain_id,
                        uniprot_id, motif_start, motif_end,
                    )
                except Exception as e:
                    logging.warning(f"UniProt enrichment failed for {apfid}: {e}")
                protein_weight = calc_protein_weight(up_data.get('protein_weight_seq', ''))
                protein_name = up_data.get('protein_name', row.get('protein_name', ''))
                biological_process = up_data.get('biological_process', row.get('biological_process', ''))
                molecular_function = up_data.get('molecular_function', row.get('molecular_function', ''))
                organism = up_data.get('organism', row.get('organism', ''))
                organism_id = up_data.get('organism_id', row.get('organism_id', ''))
                protein_partners = up_data.get('protein_partners', '')
                has_binding_sites = up_data.get('has_binding_sites', 'No')
                binding_sites = up_data.get('binding_sites', '')
            else:
                # Fallback — оставляем исходные данные
                protein_name = row.get('protein_name', '')
                biological_process = row.get('biological_process', '')
                molecular_function = row.get('molecular_function', '')
                organism = row.get('organism', '')
                organism_id = row.get('organism_id', '')
                protein_partners = ''
                has_binding_sites = 'No'
                binding_sites = ''
                protein_weight = 0.0

            # 7. weight_percent
            weight_percent = (motif_weight / protein_weight * 100) if protein_weight > 0 else 0.0
            # motif_len — из CSV или пересчитаем
            motif_len = motif_end - motif_start

            # Собираем результат
            result_dict = {
                "apfid": apfid,
                "chain_len": chain_len,
                "motif_len": motif_len,
                "fasta": fasta,
                "organism": organism,
                "organism_id": organism_id,
                "protein_name": protein_name,
                "biological_process": biological_process,
                "molecular_function": molecular_function,
                "rg": round(rg, 6),
                "sasa": round(sasa_isolated, 4),
                "surface_k": round(surface_k, 6),
                **interactions,
                "motif_weight": round(motif_weight, 4),
                "protein_weight": round(protein_weight, 4),
                "weight_percent": round(weight_percent, 6),
                "protein_partners": protein_partners,
                "has_binding_sites": has_binding_sites,
                "binding_sites": binding_sites,
                **ss_stats,
                "gravy": round(gravy, 6),
            }

            save_result(result_dict, global_idx, apfid)
            logging.info(f"OK: {apfid}")

        except Exception as e:
            err = f"Processing error: {e}"
            save_incorrect_result(apfid, err, global_idx)
            if cfg.STOP_ON_ERROR:
                raise
            log_error(apfid, sys.exc_info())
        finally:
            # Cleanup temp files
            for path in (motif_pdb_path, chain_pdb_path):
                try:
                    if os.path.exists(path):
                        os.remove(path)
                except OSError:
                    pass


# ─── Точка входа ──────────────────────────────────────────────────────────


def main():
    setup_logging()
    logging.info("Starting dbenrich26 processor")
    logging.info(
        "Config: START_INDEX=%s, STOP_INDEX=%s, ONLY_NOT_PROCESSED=%s, "
        "STOP_ON_ERROR=%s, NUM_WORKERS=%s",
        cfg.START_INDEX, cfg.STOP_INDEX, cfg.ONLY_NOT_PROCESSED,
        cfg.STOP_ON_ERROR, cfg.NUM_WORKERS,
    )

    # Инициализация
    os.makedirs(cfg.RESULTS_DIR, exist_ok=True)
    os.makedirs(cfg.TEMP_DIR, exist_ok=True)
    init_uniprot(cfg.UNIPROT_CACHE_DIR)

    # Потоковый генератор задач
    tasks = iter_experiment_groups(cfg.DATA_FILE, cfg.START_INDEX, cfg.STOP_INDEX)

    num_workers = cfg.NUM_WORKERS or os.cpu_count()
    logging.info("Starting pool with %d workers", num_workers)

    with multiprocessing.Pool(processes=num_workers) as pool:
        for _ in pool.imap(process_experiment_group, tasks):
            pass

    logging.info("Processing complete")


if __name__ == "__main__":
    main()
