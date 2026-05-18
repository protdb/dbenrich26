"""Module for calculating inter-residue interactions using ProDy."""
import logging
import os
import re
import subprocess
import tempfile
from typing import Set, Dict, Optional

import prody

import config as cfg

# Подавляем избыточный вывод ProDy в stdout
prody.confProDy(verbosity='none')


def add_hydrogens_obabel(pdb_path: str) -> Optional[str]:
    """
    Добавляет атомы водорода в PDB-файл через obabel (системный бинарник).

    Создаёт временный файл рядом с оригиналом (с суффиксом _addH.pdb).
    Возвращает путь к новому файлу или None, если obabel недоступен / завершился с ошибкой.

    obabel вызывается как subprocess — не требует Python-биндингов.
    """
    out_path = pdb_path.replace('.pdb', '_addH.pdb')
    if os.path.exists(out_path):
        return out_path
    try:
        result = subprocess.run(
            [cfg.OBABEL_BIN, pdb_path, '-O', out_path, '-h', '--quiet'],
            capture_output=True,
            timeout=30,
        )
        if result.returncode != 0 or not os.path.exists(out_path):
            logging.warning(
                "obabel failed for %s: %s",
                pdb_path,
                result.stderr.decode(errors='replace').strip(),
            )
            return None
        return out_path
    except FileNotFoundError:
        logging.warning("obabel binary not found at '%s', skipping hydrogen addition", cfg.OBABEL_BIN)
        return None
    except subprocess.TimeoutExpired:
        logging.warning("obabel timed out for %s", pdb_path)
        return None
    except Exception as e:
        logging.warning("obabel error for %s: %s", pdb_path, e)
        return None


def _parse_resnum(resname_resnum: str) -> int:
    """
    Извлекает sequence number из строки вида 'ALA42' или 'CYS123B'.
    Возвращает int, либо -1 если не удалось распарсить.
    """
    m = re.search(r'(\d+)', resname_resnum)
    return int(m.group(1)) if m else -1


def calc_interactions(chain_pdb_path: str, motif_residue_ids: Set[int]) -> Dict[str, int]:
    """
    Вычисляет межостаточные взаимодействия в цепи через ProDy.
    Разбивает на inner (оба остатка в мотиве) и outer (один в мотиве, второй снаружи).

    Если cfg.ADD_HYDROGENS=True — сначала добавляет атомы водорода через obabel
    (требуется для корректного расчёта водородных связей). При недоступности
    obabel расчёт продолжается на структуре без водородов (HB будут неполными).

    Параметры:
        chain_pdb_path: путь к PDB-файлу всей цепи
        motif_residue_ids: set sequence numbers остатков, входящих в мотив

    Возвращает словарь:
        inner_hb, inner_hp, inner_ds, outer_hb, outer_hp, outer_ds

    Формат ProDy (v2.x):
        calcHydrogenBonds / calcHydrophobic / calcDisulfideBonds возвращают
        список plain list-ов: [res1_str, atom1, chain1, res2_str, atom2, chain2, dist, ...]
        где res1_str имеет вид 'ALA42'. Метод getResnums() у этих объектов отсутствует.
    """
    # Путь к PDB, который будем подавать в ProDy
    pdb_to_use = chain_pdb_path
    h_file = None  # временный файл с водородами (удалим в finally)

    try:
        if cfg.ADD_HYDROGENS:
            h_file = add_hydrogens_obabel(chain_pdb_path)
            if h_file is not None:
                pdb_to_use = h_file

        # parsePDB возвращает AtomGroup; calcHydrogenBonds/calcHydrophobic/calcDisulfideBonds
        # в ProDy 2.x требуют Selection (у неё есть .getAtomGroup()), а не AtomGroup.
        ag = prody.parsePDB(pdb_to_use)
        structure = ag.select('all')

        result = {
            "inner_hb": 0, "inner_hp": 0, "inner_ds": 0,
            "outer_hb": 0, "outer_hp": 0, "outer_ds": 0,
        }

        # ── Hydrogen bonds ────────────────────────────────────────────────
        try:
            hbonds = prody.calcHydrogenBonds(structure)
        except Exception as e:
            logging.debug("calcHydrogenBonds failed: %s", e)
            hbonds = []
        for bond in (hbonds or []):
            # bond = ['ALA42', 'N_5', 'A', 'GLY55', 'O_12', 'A', 3.2, 25.4]
            try:
                r1 = _parse_resnum(bond[0])
                r2 = _parse_resnum(bond[3])
            except (IndexError, TypeError):
                continue
            inner1 = r1 in motif_residue_ids
            inner2 = r2 in motif_residue_ids
            if inner1 and inner2:
                result["inner_hb"] += 1
            elif inner1 or inner2:
                result["outer_hb"] += 1

        # ── Hydrophobic contacts ──────────────────────────────────────────
        try:
            hp = prody.calcHydrophobic(structure)
        except Exception as e:
            logging.debug("calcHydrophobic failed: %s", e)
            hp = []
        for contact in (hp or []):
            try:
                r1 = _parse_resnum(contact[0])
                r2 = _parse_resnum(contact[3])
            except (IndexError, TypeError):
                continue
            inner1 = r1 in motif_residue_ids
            inner2 = r2 in motif_residue_ids
            if inner1 and inner2:
                result["inner_hp"] += 1
            elif inner1 or inner2:
                result["outer_hp"] += 1

        # ── Disulfide bonds ───────────────────────────────────────────────
        try:
            ds = prody.calcDisulfideBonds(structure)
        except Exception as e:
            logging.debug("calcDisulfideBonds failed: %s", e)
            ds = []
        for bond in (ds or []):
            try:
                r1 = _parse_resnum(bond[0])
                r2 = _parse_resnum(bond[3])
            except (IndexError, TypeError):
                continue
            inner1 = r1 in motif_residue_ids
            inner2 = r2 in motif_residue_ids
            if inner1 and inner2:
                result["inner_ds"] += 1
            elif inner1 or inner2:
                result["outer_ds"] += 1

        return result

    finally:
        # Удаляем временный файл с водородами
        if h_file and os.path.exists(h_file):
            try:
                os.remove(h_file)
            except OSError:
                pass
