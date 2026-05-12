import json
from gzip import GzipFile
from time import sleep

import jmespath
import pandas as pd
import os
import csv
from uniprot_meta_tool import UniprotData
from uniprot_meta_tool.data_parser import set_config
from get_pdb_api_data import get_single
from apfid import Apfid
from mass_calc import calc_mass
from pathlib import Path

import logging

logging.getLogger().setLevel(logging.DEBUG)

set_config(meta_cache_dir='/home/kirill/meta/uniprot/')

out_path = 'metadata.csv'
file_exists = os.path.exists(out_path)
if file_exists:
    with open(out_path, 'r') as f:
        reader = csv.DictReader(f)
        done_apfids = [x['apfid'] for x in reader]
else:
    done_apfids = []

# apfids = [x for x in pd.read_csv('in_structs/data_all.csv')['apfid'] if x not in done_apfids]

DATA_PLAIN = "data.entry.polymer_entities[*].{chains: polymer_entity_instances[*].rcsb_polymer_entity_instance_container_identifiers.auth_asym_id, uniprot_id: uniprots[0].rcsb_id, protein_name:uniprots[0].rcsb_uniprot_protein.name.value, organism_id: rcsb_entity_source_organism[0].ncbi_taxonomy_id,organism:rcsb_entity_source_organism[0].ncbi_scientific_name, type:entity_poly.rcsb_entity_polymer_type}"

metadata_dir = Path('/home/kirill/meta/pdb_single')

with open(out_path, 'a') as out, open('in_structs/data_all.csv') as in_file:
    writer = None
    reader = csv.DictReader(in_file)
    for line in reader:
        struct_id = line['apfid']
        if struct_id in done_apfids:
            continue
        apfid = Apfid(apfid=struct_id)
        file_path = metadata_dir / f'{apfid.experiment_id}.json.gz'
        if file_path.exists():
            if file_path.stat().st_size > 0:
                with GzipFile(file_path, 'rb') as f:
                    plain_data = json.load(f)
        else:
            try:
                apfid = Apfid(apfid=struct_id)
                print(apfid.upper())
                data = get_single(apfid.experiment_id)
            except Exception as e:
                logging.exception(e)
                sleep(1)
                continue
            try:
                plain_data = [x for x in jmespath.search(DATA_PLAIN, data) if x['type'] == 'Protein']
            except TypeError:
                continue
        try:
            up_id = [x['uniprot_id'] for x in plain_data if apfid.chain_id in x['chains']][0]
            print(up_id)
        except IndexError:
            continue
        partner_uniprots = list(set([x['uniprot_id'] for x in plain_data if x['uniprot_id'] != up_id]))
        if len(partner_uniprots) > 0:
            try:
                partner_ups = [UniprotData(x) for x in partner_uniprots if type(x) is str]
                protein_partners = [f"{x.gene} ({x.uniprot_id})" if x.gene else f"{x.uniprot_id}" for x in partner_ups]
            except TypeError:
                protein_partners = []
        else:
            protein_partners = []
        protein_partners = '; '.join(protein_partners)
        # print(protein_partners)
        if str(up_id) == 'None':
            continue
        try:
            up_meta = UniprotData(up_id)
        except:
            continue
        molecular_weight = calc_mass(up_meta.sequence)
        if apfid.start is not None:
            try:
                features = up_meta.get_position_features(
                    apfid.start,
                    apfid.end,
                )
            except TypeError:
                features = []
        else:
            features = up_meta.features
        binding_site_features = [x for x in features if x.feature_type == 'Binding site']
        has_binding_sites = len(binding_site_features) > 0
        if has_binding_sites:
            binding_sites = '; '.join(
                [str(x.start) if x.start == x.end else f'{x.start}-{x.end}' for x in binding_site_features]
            )
        else:
            binding_sites = ''
        result = {
            'apfid': struct_id,
            'protein_weight': molecular_weight,
            'protein_partners': protein_partners,
            'has_binding_sites': has_binding_sites,
            'binding_sites': binding_sites,
        }
        if writer is None:
            writer = csv.DictWriter(out, fieldnames=result.keys())
            if not file_exists:
                writer.writeheader()
        writer.writerow(result)
        out.flush()