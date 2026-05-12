"""Module for calculating inter-residue interactions using ProDy."""
import prody
from typing import Set, Dict


def calc_interactions(chain_pdb_path: str, motif_residue_ids: Set[int]) -> Dict[str, int]:
    """
    Вычисляет межостаточные взаимодействия в цепи через ProDy.
    Разбивает на inner (оба остатка в мотиве) и outer (один в мотиве, второй снаружи).

    Параметры:
        chain_pdb_path: путь к PDB-файлу всей цепи
        motif_residue_ids: set sequence numbers остатков, входящих в мотив

    Возвращает словарь:
        inner_hb, inner_hp, inner_ds, outer_hb, outer_hp, outer_ds
    """
    # Загружаем цепь
    structure = prody.parsePDB(chain_pdb_path)
    inter = prody.Interactions()

    hbonds = inter.calcHydrogenBonds(structure)
    hp = inter.calcHydrophobicContacts(structure)
    ds = inter.calcDisulfideBonds(structure)

    result = {
        "inner_hb": 0, "inner_hp": 0, "inner_ds": 0,
        "outer_hb": 0, "outer_hp": 0, "outer_ds": 0,
    }

    # Hydrogen bonds
    if hbonds is not None:
        for bond in hbonds:
            # bond — это prody объект связи, содержит resnum участников
            try:
                r1, r2 = bond.getResnums()
            except (AttributeError, ValueError):
                continue
            inner1 = r1 in motif_residue_ids
            inner2 = r2 in motif_residue_ids
            if inner1 and inner2:
                result["inner_hb"] += 1
            elif inner1 or inner2:
                result["outer_hb"] += 1

    # Hydrophobic contacts
    if hp is not None:
        for contact in hp:
            try:
                r1, r2 = contact.getResnums()
            except (AttributeError, ValueError):
                continue
            inner1 = r1 in motif_residue_ids
            inner2 = r2 in motif_residue_ids
            if inner1 and inner2:
                result["inner_hp"] += 1
            elif inner1 or inner2:
                result["outer_hp"] += 1

    # Disulfide bonds
    if ds is not None:
        for bond in ds:
            try:
                r1, r2 = bond.getResnums()
            except (AttributeError, ValueError):
                continue
            inner1 = r1 in motif_residue_ids
            inner2 = r2 in motif_residue_ids
            if inner1 and inner2:
                result["inner_ds"] += 1
            elif inner1 or inner2:
                result["outer_ds"] += 1

    return result
