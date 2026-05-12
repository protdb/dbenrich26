"""Module for SASA calculations using FreeSASA with in-memory BioPython structures."""
import freesasa
from Bio.PDB import Structure, Model, Chain, Residue
from typing import List


def _build_single_chain_structure(residues: List[Residue], structure_id: str = "tmp") -> Structure.Structure:
    """
    Строит временную структуру из списка остатков (все в одну цепь 'A')
    для передачи в FreeSASA. Используется только в памяти, не сохраняется на диск.
    """
    struct = Structure.Structure(structure_id)
    model = Model.Model(0)
    chain = Chain.Chain('A')
    for res in residues:
        # Копируем остаток с его атомами
        new_res = Residue.Residue(res.id, res.resname, res.segid)
        for atom in res:
            new_res.add(atom.copy())
        chain.add(new_res)
    model.add(chain)
    struct.add(model)
    return struct


def calc_sasa_isolated(motif_residues: List[Residue]) -> float:
    """
    Вычисляет SASA мотива в изоляции (только остатки мотива, без контекста цепи).
    Использует структуру в памяти — без записи на диск.
    Возвращает суммарную SASA в Å².
    """
    tmp_struct = _build_single_chain_structure(motif_residues)
    result = freesasa.calcBioPDB(tmp_struct)
    return result.totalArea()


def calc_sasa_in_chain(chain, motif_residue_ids: List[int]) -> float:
    """
    Вычисляет SASA для полной цепи, но суммирует только по остаткам,
    входящим в мотив (по sequence number из motif_residue_ids).
    Возвращает SASA этих остатков в контексте цепи.
    """
    # Собираем структуру из всей цепи
    all_residues = list(chain)
    tmp_struct = _build_single_chain_structure(all_residues)
    result = freesasa.calcBioPDB(tmp_struct)

    # Суммируем SASA только по атомам нужных остатков
    target = set(motif_residue_ids)
    total = 0.0
    for atom_result in result.atoms():
        # atom_result имеет атрибуты: residueNum, area
        if atom_result.residueNum in target:
            total += atom_result.area
    return total


def calc_surface_k(sasa_isolated: float, sasa_in_chain: float) -> float:
    """
    Коэффициент погружённости мотива:
    surface_k = SASA_в_цепи / SASA_в_изоляции
    Если SASA_в_изоляции == 0, возвращает 0.0.
    """
    if sasa_isolated == 0.0:
        return 0.0
    return sasa_in_chain / sasa_isolated
