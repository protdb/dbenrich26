"""Module for UniProt metadata retrieval and enrichment."""
import json
from gzip import GzipFile
from pathlib import Path
from typing import Optional, Dict

import jmespath
from uniprot_meta_tool import UniprotData
from uniprot_meta_tool.data_parser import set_config
from uniprot_meta_tool.utils import get_biological_processes, get_molecular_functions
from get_pdb_api_data import get_single
from apfid import Apfid

# JMESPath запрос для получения данных PDB entry (из примера)
DATA_PLAIN = (
    "data.entry.polymer_entities[*].{chains: polymer_entity_instances[*]."
    "rcsb_polymer_entity_instance_container_identifiers.auth_asym_id, "
    "uniprot_id: uniprots[0].rcsb_id, "
    "protein_name:uniprots[0].rcsb_uniprot_protein.name.value, "
    "organism_id: rcsb_entity_source_organism[0].ncbi_taxonomy_id, "
    "organism:rcsb_entity_source_organism[0].ncbi_scientific_name, "
    "type:entity_poly.rcsb_entity_polymer_type}"
)


def init_uniprot(cache_dir: str) -> None:
    """Инициализирует конфигурацию UniProt (кэш-директорию)."""
    set_config(meta_cache_dir=cache_dir)


def get_uniprot_data(uniprot_id: str) -> Optional[Dict]:
    """
    Получает метаданные белка из UniProt по ID.
    Возвращает словарь с полями:
        protein_name, biological_process, molecular_function,
        protein_weight (sequence string — для последующего расчёта массы),
        organism, organism_id, protein_partners, has_binding_sites, binding_sites

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
        'protein_partners': '',
        'has_binding_sites': 'No',
        'binding_sites': '',
    }


def enrich_with_partners_and_sites(
    up_data: Dict,
    experiment_id: str,
    chain_id: str,
    uniprot_id: str,
    motif_start: int,
    motif_end: int,
) -> Dict:
    """
    Дополняет словарь up_data полями protein_partners, has_binding_sites, binding_sites.
    Использует PDB API для получения партнёров и UniprotData.get_position_features
    для сайтов связывания.

    Аргументы:
        up_data: результат get_uniprot_data (будет изменён)
        experiment_id: PDB ID
        chain_id: идентификатор цепи
        uniprot_id: UniProt ID
        motif_start, motif_end: координаты мотива

    Возвращает тот же словарь (up_data) с заполненными полями.
    """
    # Партнёры по белку — через PDB API
    try:
        apfid_obj = Apfid(experiment_id=experiment_id, chain_id=chain_id)
    except Exception:
        apfid_obj = None

    plain_data = None
    if apfid_obj:
        try:
            data = get_single(experiment_id)
            plain_data = [
                x for x in jmespath.search(DATA_PLAIN, data)
                if x and x.get('type') == 'Protein'
            ]
        except Exception:
            plain_data = None

    if plain_data:
        try:
            # Ищем запись с нашей цепью
            current_entry = None
            for entry in plain_data:
                chains = entry.get('chains') or []
                if chain_id in chains:
                    current_entry = entry
                    break
            if current_entry:
                partner_uniprots = list(set([
                    x['uniprot_id'] for x in plain_data
                    if x['uniprot_id'] != uniprot_id and x.get('uniprot_id')
                ]))
                if partner_uniprots:
                    try:
                        partner_ups = [UniprotData(uid) for uid in partner_uniprots if isinstance(uid, str)]
                        protein_partners_list = [
                            f"{p.gene} ({p.uniprot_id})" if p.gene else p.uniprot_id
                            for p in partner_ups
                        ]
                        up_data['protein_partners'] = '; '.join(protein_partners_list)
                    except Exception:
                        up_data['protein_partners'] = ''
                else:
                    up_data['protein_partners'] = ''
        except Exception:
            up_data['protein_partners'] = ''

    # Сайты связывания
    try:
        up = UniprotData(uniprot_id)
    except Exception:
        up_data['has_binding_sites'] = 'No'
        up_data['binding_sites'] = ''
        return up_data

    try:
        features = up.get_position_features(motif_start, motif_end)
    except Exception:
        features = []

    binding_site_features = [f for f in features if f.feature_type == 'Binding site']
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
