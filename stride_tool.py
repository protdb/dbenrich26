"""Module for running STRIDE and parsing secondary structure assignments."""
import subprocess
import os
from typing import Optional, List, Dict

from config import COIL_TAIL_MAX_LEN

STRIDE_BIN = "./stride"


def run_stride(pdb_path: str) -> Optional[str]:
    """
    Запускает STRIDE на PDB-файле.
    Возвращает raw-текст результата.
    Временный .stride-файл удаляется после чтения.
    При ошибке возвращает None.
    """
    inp = pdb_path
    outp = pdb_path + '.stride'
    try:
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(
                [STRIDE_BIN, '-f' + outp, inp],
                stderr=devnull, stdout=devnull
            )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    try:
        with open(outp) as f:
            raw = f.read()
    except FileNotFoundError:
        return None
    finally:
        if os.path.exists(outp):
            os.remove(outp)
    return raw


def parse_stride_output(raw: str) -> Optional[List[Dict]]:
    """
    Парсит ASG-строки из вывода STRIDE.
    Возвращает список словарей с ключами:
        amino, total_no, chain, sec_letter, sec_name
    Ремаппинг: T, B, G → C (coil).
    """
    if not raw:
        return None
    try:
        lines = raw.splitlines(keepends=False)
        asg_lines = [l for l in lines if l[:3] == "ASG"]
    except AttributeError:
        return None

    parsed = []
    for line in asg_lines:
        try:
            sec_letter = line[24]
            if sec_letter in ('T', 'B', 'G'):
                sec_letter = 'C'
            parsed.append({
                "amino": line[5:8].strip(),
                "total_no": int(line[10:15]),
                "chain": line[9].strip(),
                "sec_letter": sec_letter,
                "sec_name": line[25:39].strip(),
            })
        except (ValueError, IndexError):
            continue
    return parsed if parsed else None


def get_secondary_structure_elements(parsed: List[Dict]) -> List[Dict]:
    """
    Схлопывает последовательные остатки с одинаковым sec_letter
    в элементы вторичной структуры.
    Каждый элемент: {type, sec_letter, start, end, length}
    """
    if not parsed:
        return []
    elements = []
    current = {
        "type": parsed[0]["sec_name"],
        "sec_letter": parsed[0]["sec_letter"],
        "start": parsed[0]["total_no"],
    }
    for i in range(1, len(parsed)):
        if parsed[i]["sec_letter"] != current["sec_letter"]:
            current["end"] = parsed[i - 1]["total_no"]
            current["length"] = current["end"] - current["start"] + 1
            elements.append(current)
            current = {
                "type": parsed[i]["sec_name"],
                "sec_letter": parsed[i]["sec_letter"],
                "start": parsed[i]["total_no"],
            }
    current["end"] = parsed[-1]["total_no"]
    current["length"] = current["end"] - current["start"] + 1
    elements.append(current)
    return elements


def get_motif_ss_stats(
    elements: List[Dict],
    motif_start: int,
    motif_end: int,
    coil_tail_max_len: int = COIL_TAIL_MAX_LEN,
) -> Dict[str, int]:
    """
    Фильтрует SS-элементы, полностью входящие в [motif_start, motif_end].
    Удаляет хвостовые coil'ы (первый/последний не длиннее coil_tail_max_len).
    Возвращает статистику: helix/strand/coil count + длины.
    """
    # Фильтруем: элемент должен полностью лежать внутри мотива
    filtered = [
        e for e in elements
        if e["start"] >= motif_start and e["end"] <= motif_end
    ]

    # Удаляем хвостовые coil'ы (однократно, не рекурсивно)
    if filtered and filtered[0]["sec_letter"] == 'C' and filtered[0]["length"] <= coil_tail_max_len:
        filtered = filtered[1:]
    if filtered and filtered[-1]["sec_letter"] == 'C' and filtered[-1]["length"] <= coil_tail_max_len:
        filtered = filtered[:-1]

    helix_count = 0
    strand_count = 0
    coil_count = 0
    helix_length = 0
    strand_length = 0
    coil_length = 0

    for e in filtered:
        if e["sec_letter"] in ('H', 'I', 'G'):
            helix_count += 1
            helix_length += e["length"]
        elif e["sec_letter"] == 'E':
            strand_count += 1
            strand_length += e["length"]
        elif e["sec_letter"] == 'C':
            coil_count += 1
            coil_length += e["length"]

    return {
        "helix_count": helix_count,
        "strand_count": strand_count,
        "coil_count": coil_count,
        "helix_length": helix_length,
        "strand_length": strand_length,
        "coil_length": coil_length,
    }
