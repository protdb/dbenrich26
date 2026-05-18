"""Module for UniProt metadata retrieval and enrichment."""
from typing import Optional, Dict

from uniprot_meta_tool import UniprotData
from uniprot_meta_tool.data_parser import set_config
from uniprot_meta_tool.utils import get_biological_processes, get_molecular_functions


def init_uniprot(cache_dir: str) -> None:
    """Инициализирует конфигурацию UniProt (кэш-директорию)."""
    set_config(meta_cache_dir=cache_dir)


def get_uniprot_data(uniprot_id: str) -> Optional[Dict]:
    """
    Получает метаданные белка из UniProt по ID.
    Возвращает словарь с полями:
        protein_name, biological_process, molecular_function,
        protein_weight_seq (sequence string — для последующего расчёта массы),
        organism, organism_id

    При ошибке (пустой ID, сетевая ошибка, obsolete) возвращает None.
    """
    if not uniprot_id or uniprot_id.strip() == '':
        return None

    try:
        up = UniprotData(uniprot_id)
    except Exception:
        return None

    return {
        'protein_name': up.name or '',
        'biological_process': '; '.join(get_biological_processes(up)) if up.keywords else '',
        'molecular_function': '; '.join(get_molecular_functions(up)) if up.keywords else '',
        'protein_weight_seq': up.sequence or '',
        'organism': up.organism.scientific_name if up.organism else '',
        'organism_id': up.organism.taxon_id if up.organism else '',
    }


def enrich_binding_sites(
    up_data: Dict,
    uniprot_id: str,
    motif_start: int,
    motif_end: int,
) -> Dict:
    """
    Дополняет словарь up_data полями protein_partners, has_binding_sites, binding_sites
    на основе UniProt-аннотации (get_position_features).

    Поле protein_partners всегда пустое — получение партнёров через PDB API
    в текущей версии не используется.

    Аргументы:
        up_data: результат get_uniprot_data (будет изменён)
        uniprot_id: UniProt ID
        motif_start, motif_end: координаты мотива

    Возвращает тот же словарь (up_data) с заполненными полями.
    """
    # Сайты связывания — через UniProt
    try:
        up = UniprotData(uniprot_id)
    except Exception:
        up_data['protein_partners'] = ''
        up_data['has_binding_sites'] = 'No'
        up_data['binding_sites'] = ''
        return up_data

    try:
        features = up.get_position_features(motif_start, motif_end)
    except Exception:
        features = []

    binding_site_features = [f for f in features if f.feature_type == 'Binding site']
    up_data['protein_partners'] = ''
    up_data['has_binding_sites'] = 'Yes' if binding_site_features else 'No'
    if binding_site_features:
        parts = []
        for bs in binding_site_features:
            if bs.start == bs.end:
                parts.append(str(bs.start))
            else:
                parts.append(f'{bs.start}-{bs.end}')
        up_data['binding_sites'] = '; '.join(parts)
    else:
        up_data['binding_sites'] = ''

    return up_data