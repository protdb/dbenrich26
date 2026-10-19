from python_graphql_client.graphql_client import GraphqlClient


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


def get_single(pdb_id):
    gql_client = GraphqlClient("https://data.rcsb.org/graphql")
    return gql_client.execute(
        query=PDB_TEMPLATE,
        variables={
            'pdb_id': pdb_id
        }
    )