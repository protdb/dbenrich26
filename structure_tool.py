"""Module for PDB structure loading, residue manipulation, and basic calculations."""
import gzip
import math
import os
from typing import List, Optional

from Bio.PDB import PDBParser, PDBIO, Select, Structure
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.SeqUtils import molecular_weight

from config import PDB_DIR, TEMP_DIR, MISSING_RESIDUE_THRESHOLD


# ─── Шкала Kyte-Doolittle для GRAVY ─────────────────────────────────────────
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
}
AVG_RESIDUE_MASS = 111.1  # средняя масса аминокислотного остатка (Да)
MASS_OF_ALANINE = 89.093  # масса Ala для коррекции при замене X → A


class MotifValidationError(Exception):
    """Бросается, когда мотив не проходит проверку (пустой / слишком много пропусков)."""
    pass


# ─── PDB путь и загрузка ────────────────────────────────────────────────────


def get_pdb_path(experiment_id: str) -> str:
    """
    Конструирует путь к PDB-файлу в архиве.
    Пример: '4E06' → '{PDB_DIR}/e0/pdb4e06.ent.gz'
    """
    eid = experiment_id.lower()
    middle = eid[1:3]
    return os.path.join(PDB_DIR, middle, f"pdb{eid}.ent.gz")


def load_structure(experiment_id: str) -> Structure.Structure:
    """
    Загружает PDB-структуру из gz-архива через BioPython.
    Возвращает объект Structure (первая модель).
    Бросает FileNotFoundError если файл не найден.
    """
    path = get_pdb_path(experiment_id)
    if not os.path.exists(path):
        raise FileNotFoundError(f"PDB file not found: {path}")
    parser = PDBParser(QUIET=True)
    with gzip.open(path, 'rt') as f:
        structure = parser.get_structure(experiment_id, f)
    return structure


# ─── Фильтрация остатков ────────────────────────────────────────────────────


def is_protein_residue(residue) -> bool:
    """
    Проверяет, является ли остаток белковым (стандартная или модифицированная аминокислота).
    - Стандартные: hetflag == ' '
    - Модифицированные в основной цепи: hetflag == ' ' или hetflag.startswith('H_') и есть CA-атом
    - Исключаем: воду (HOH), лиганды, ионы
    """
    if residue.get_resname().strip() == 'HOH':
        return False
    hetflag = residue.id[0]
    if hetflag == ' ':
        return True
    if hetflag.startswith('H_') and 'CA' in residue:
        return True
    return False


def get_protein_residues(chain) -> list:
    """Возвращает список белковых остатков цепи (только residue с is_protein_residue)."""
    return [res for res in chain if is_protein_residue(res)]


def get_chain_len(chain) -> int:
    """Возвращает количество белковых остатков в цепи."""
    return len(get_protein_residues(chain))


# ─── Вырезание мотива ──────────────────────────────────────────────────────


def extract_motif_residues(chain, motif_start: int, motif_end: int) -> list:
    """
    Из белковых остатков цепи берёт те, чей sequence number (residue.id[1])
    в диапазоне [motif_start, motif_end] включительно.
    iCode (residue.id[2]) игнорируется.
    """
    protein_residues = get_protein_residues(chain)
    return [
        res for res in protein_residues
        if motif_start <= res.id[1] <= motif_end
    ]


def validate_motif(motif_residues: list, motif_start: int, motif_end: int) -> None:
    """
    Проверяет мотив:
    - Не пустой
    - Количество остатков >= ожидаемое * (1 - MISSING_RESIDUE_THRESHOLD)
    """
    if not motif_residues:
        raise MotifValidationError("Motif is empty (no residues found in range)")
    expected = motif_end - motif_start
    actual = len(motif_residues)
    if actual < expected * (1 - MISSING_RESIDUE_THRESHOLD):
        raise MotifValidationError(
            f"Too many missing residues: expected ~{expected}, got {actual} "
            f"(threshold: {MISSING_RESIDUE_THRESHOLD*100:.0f}%)"
        )


# ─── Сохранение PDB-фрагментов (для STRIDE / ProDy) ────────────────────────


class _ResidueSelect(Select):
    """Bio.PDB Select class для фильтрации по заданному набору остатков."""
    def __init__(self, residue_set: set):
        self.residue_set = residue_set

    def accept_residue(self, residue) -> bool:
        return residue in self.residue_set


def save_motif_to_pdb(motif_residues: list, chain_id: str, out_path: str) -> None:
    """
    Сохраняет остатки мотива во временный PDB-файл.
    Используется для запуска STRIDE.
    """
    # Берём родительскую структуру первого остатка
    chain = motif_residues[0].get_parent()
    model = chain.get_parent()
    structure = model.get_parent()
    residue_set = set(motif_residues)
    io = PDBIO()
    io.set_structure(structure)
    io.save(out_path, select=_ResidueSelect(residue_set))


def save_chain_to_pdb(chain, out_path: str) -> None:
    """
    Сохраняет всю цепь во временный PDB-файл.
    Используется для ProDy.
    """
    io = PDBIO()
    io.set_structure(chain.get_parent().get_parent())
    io.save(out_path, select=_ResidueSelect(set(get_protein_residues(chain))))


# ─── Радиус гирации ─────────────────────────────────────────────────────────


def calc_rg(motif_residues: list) -> float:
    """
    Радиус гирации по Cα-атомам мотива.
    Rg = sqrt( sum((r_i - r_cm)^2) / N )
    """
    coords = []
    for res in motif_residues:
        if 'CA' in res:
            coords.append(res['CA'].get_coord())
    if not coords:
        raise MotifValidationError("No CA atoms found in motif for Rg calculation")
    center = sum(coords) / len(coords)
    sq_dist = sum(((c - center) ** 2).sum() for c in coords)
    return math.sqrt(sq_dist / len(coords))


# ─── Молекулярная масса ─────────────────────────────────────────────────────


def get_motif_fasta(motif_residues: list) -> str:
    """
    Строит однобуквенную последовательность мотива.
    Нестандартные остатки → '?'.
    """
    chars = []
    for res in motif_residues:
        try:
            chars.append(protein_letters_3to1[res.resname])
        except KeyError:
            chars.append('?')
    return ''.join(chars)


def calc_motif_weight(motif_residues: list) -> float:
    """
    Молекулярная масса мотива через Bio.SeqUtils.molecular_weight.
    '?' заменяются на 'A' с коррекцией на среднюю массу.
    """
    fasta = get_motif_fasta(motif_residues)
    return calc_protein_weight(fasta)


def calc_protein_weight(sequence: str) -> float:
    """
    Молекулярная масса последовательности.
    '?' / 'X' → средняя масса остатка 111.1 Да.
    Если последовательность пуста → 0.0.
    """
    if not sequence:
        return 0.0
    # Заменяем X/? на A, BioPython их не понимает
    clean_seq = sequence.replace('?', 'A').replace('X', 'A')
    try:
        weight = molecular_weight(clean_seq, seq_type='protein')
    except Exception:
        return 0.0
    # Коррекция: 'A' весит ~89.09, а мы хотим 111.1 для каждого заменённого
    x_count = sequence.count('?') + sequence.count('X')
    if x_count > 0:
        weight += x_count * (AVG_RESIDUE_MASS - MASS_OF_ALANINE)
    return weight


# ─── GRAVY ──────────────────────────────────────────────────────────────────


def calc_gravy(fasta: str) -> float:
    """
    GRAVY (Grand Average of Hydropathy) — средняя гидрофобность
    по шкале Kyte-Doolittle.
    '?' и 'X' пропускаются (не участвуют в среднем).
    Если нет ни одного известного остатка → 0.0.
    """
    scores = [KD_SCALE[aa] for aa in fasta if aa in KD_SCALE]
    return sum(scores) / len(scores) if scores else 0.0
