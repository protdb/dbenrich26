import logging
from functools import lru_cache
from apfid import Apfid
from python_graphql_client.graphql_client import GraphqlClient
import jmespath
from uniprot_meta_tool import UniprotData

UNIPROT_TEMPLATE = """
query ($pdb_ids: [String!]!) {
  entries(entry_ids: $pdb_ids) {
    rcsb_id
    polymer_entities {
      rcsb_id
      polymer_entity_instances {
        rcsb_id
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
          entry_id
        }
      }
      uniprots {
        rcsb_id
      }
    }
  }
}
"""

PDB_TEMPLATE = """
query structure ($pdb_id: String!) {
  entry(entry_id:$pdb_id){
      struct {
          title
      }
    rcsb_id
    struct_keywords {
      text
    }
    struct {
      title
    }
    rcsb_entry_container_identifiers {
      entry_id
    }
    rcsb_entry_info {
      experimental_method
      resolution_combined
    }
    pdbx_molecule_features {
      class
      details
    }

    polymer_entities {
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id
        }
      }
      rcsb_polymer_entity {
        details
      }
      uniprots {
        rcsb_id
        rcsb_uniprot_protein {
          name {
            value
          }
        }
      }
      entity_poly {
        rcsb_entity_polymer_type
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        ncbi_scientific_name
      }
    }
  }
}
"""


def get_chains(pdb_ids):
    gql_client = GraphqlClient("https://data.rcsb.org/graphql")
    return gql_client.execute(
        query=UNIPROT_TEMPLATE,
        variables={
            'pdb_ids':pdb_ids
        }
    )

@lru_cache(maxsize=16)
def get_single(pdb_id):
    gql_client = GraphqlClient("https://data.rcsb.org/graphql")
    return gql_client.execute(
        query=PDB_TEMPLATE,
        variables={
            'pdb_id': pdb_id
        }
    )


DATA_PLAIN = "data.entry.polymer_entities[*].{chains: polymer_entity_instances[*].rcsb_polymer_entity_instance_container_identifiers.auth_asym_id, uniprot_id: uniprots[0].rcsb_id, protein_name:uniprots[0].rcsb_uniprot_protein.name.value, organism_id: rcsb_entity_source_organism[0].ncbi_taxonomy_id,organism:rcsb_entity_source_organism[0].ncbi_scientific_name, type:entity_poly.rcsb_entity_polymer_type}"

params_queries = {
    "experiment_type": "data.entry.rcsb_entry_info.experimental_method",
    "resolution": "data.entry.rcsb_entry_info.resolution_combined",
    "experiment_name": "data.entry.struct.title"

}

def get_partner_uniprots(up_id, plain_data):
    partner_uniprots = list(set([x['uniprot_id'] for x in plain_data if x['uniprot_id'] != up_id]))
    if len(partner_uniprots) > 0:
        try:
            partner_ups = [UniprotData(x) for x in partner_uniprots if type(x) is str]
            protein_partners = [f"{x.gene} ({x.uniprot_id})" if x.gene else f"{x.uniprot_id}" for x in partner_ups]
        except TypeError:
            protein_partners = []
    else:
        protein_partners = []
    return '; '.join(protein_partners)

@lru_cache(maxsize=16)
def get_apfid_data(struct_id: str | Apfid)-> dict | None:
    if type(struct_id) is str:
        struct_id = Apfid(apfid=struct_id)
    print(struct_id)
    print(struct_id.experiment_id)
    try:
        data = get_single(struct_id.experiment_id)
        plain_data = jmespath.search(
            DATA_PLAIN,
            data
        )
    except Exception as e:
        logging.error(e)
        return None
    if plain_data is None:
        return None
    try:
        up_id = up_id = [x['uniprot_id'] for x in plain_data if struct_id.chain_id in x['chains']][0]
    except Exception as e:
        logging.error(e)
        up_id = ''
    all_chains = jmespath.search("[*].chains", plain_data)
    if all_chains is not None and len(all_chains) != 1:
        single_chain = False
    else:
        single_chain = len(all_chains[0]) == 1
    res = {
        'uniprot_id': up_id,
        'single_chain': single_chain,
        'protein_partners': get_partner_uniprots(up_id, plain_data),
    }
    res.update({
        k: jmespath.search(v, data) for k, v in params_queries.items()
    })
    return res
