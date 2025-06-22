import requests
import json
import pandas as pd
import re

"""
I will add:
- genes for which the lead variant is missense;
- genes for which the lead variants is in a qtl credible set;
- if none of the above,  nearest gene

# TODO query for nearest gene using variant targets?
"""
base_url = "https://api.platform.opentargets.org/api/v4/graphql"

def query_snp_id_by_rsid(rsid):
    variant_id_query = """
    query variant_id_query($rsid: String!) {
        search(queryString: $rsid, entityNames: ["variant"]) {
            hits {
                    id
                    score
                }
        }
    }
    """

    variables = {"rsid": rsid}

    r = requests.post(base_url, json = {"query": variant_id_query, "variables": variables})

    response_data = json.loads(r.text)['data']

    # return response_data

    if response_data['search']['hits'] is None or len(response_data['search']['hits']) == 0:
        return None

    daf = pd.json_normalize(response_data['search']['hits']).sort_values(by = 'score', ascending = False)

    return daf.loc[0, 'id']

def query_snp_credible_sets(variant_id, size = 1000, index = 0):
    credible_sets_query = """
    query QTLCredibleSetsQuery($variant_id: String!, $size: Int!, $index: Int!) {
    variant(variantId: $variant_id) {
        qtlCredibleSets: credibleSets(
        studyTypes: [scsqtl, sceqtl, scpqtl, sctuqtl, sqtl, eqtl, pqtl, tuqtl]
        page: { size: $size, index: $index }
        ) {
        count
        rows {
            finemappingMethod
            confidence
            isTransQtl
            variant {
            id
            }
            study {
            id
            studyType
            condition
            target {
                id
                approvedSymbol
            }
            biosample {
                biosampleId
                biosampleName
            }
            }
            locus(variantIds: [$variant_id]) {
            rows {
                posteriorProbability
                beta
                is95CredibleSet
            }
            }
            locusSize: locus {
            count
            }
        }
        }
    }
    }
    """

    variables = {"variant_id": variant_id, "size": size, "index": index}

    r = requests.post(base_url, json = {"query": credible_sets_query, "variables": variables})

    response_data = json.loads(r.text)['data']

    daf = pd.json_normalize(response_data['variant']['qtlCredibleSets']['rows'])

    # TODO handle empty daf

    daf['locus.rows_extracted'] = daf['locus.rows'].apply(lambda x: x[0] if isinstance(x, list) and x else {})

    expanded = daf["locus.rows_extracted"].apply(pd.Series)

    daf = daf.drop(["locus.rows","locus.rows_extracted"], axis=1).join(expanded)

    return daf

missense_query = """
query missense_and_credible_sets($variant_id: String!) {
    variant(variantId: $variant_id) {
    id
    rsIds
    transcriptConsequences {
        target {
        approvedSymbol
        symbolSynonyms {
            label
        }
        }
        consequenceScore
        variantConsequences {
        label
        id
        }
    }
    }
}
"""
# # NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints
# genes_query = """
# query nearest_gene_and_top_genes($variant_id: String!) {
#     variantInfo(variantId: $variant_id) {
#         nearestGene {
#             symbol
#         }
#         nearestGeneDistance
#         mostSevereConsequence
#         }
#     genesForVariant(variantId: $variant_id) {
#         overallScore
#         gene {
#         symbol
#         id
#         }
#   }
# }
# """

# base_url = "https://api.genetics.opentargets.org/graphql"

# daf = pd.read_csv(snakemake.input[0], sep = '\\s+')

# if daf.shape[0] == 0:
#     pd.DataFrame(columns = ['rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'nearest_gene', 'nearest_gene_distance', 'top_gene', 'most_severe_consequence']).to_csv(snakemake.output[0], sep = '\t')
# else:
#     res_dicts = []

#     for index, row in daf.iterrows():
#          res_dicts.append(query_snp(row.rsid, row.chromosome, row.base_pair_location, row.other_allele, row.effect_allele))

#     annot_daf = pd.DataFrame(res_dicts)

#     pd.merge(daf, annot_daf, on = 'rsid').to_csv(snakemake.output[0], sep = '\t', index = False)
