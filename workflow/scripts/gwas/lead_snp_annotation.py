import requests
import json
import pandas as pd
import re

def query_snp(rsid, chrom, bp, ref, alt):
    variables = {"variant_id": '_'.join([str(x) for x in [chrom, bp, ref, alt]])}

    r = requests.post(base_url, json={"query": genes_query, "variables": variables})

    response_data = json.loads(r.text)['data']

    if not response_data['variantInfo']:
        variables = {"variant_id": '_'.join([str(x) for x in [chrom, bp, ref, alt]])}

        r = requests.post(base_url, json={"query": genes_query, "variables": variables})

        if not response_data['variantInfo']:
            return {'rsid' : rsid,
                    'nearest_gene': '',
                    'nearest_gene_distance': '',
                    'top_genes': '',
                    'most_severe_consequence': ''
                    }

    if response_data['genesForVariant']:
        top_genes_daf = pd.json_normalize(response_data['genesForVariant']).sort_values(by = 'overallScore', ascending = False)

        top_genes_str = ','.join(top_genes_daf.head(snakemake.params.no_of_top_genes)['gene.symbol'].astype(str).tolist())
    else:
        top_genes_str = ''

    if response_data['variantInfo']['nearestGene']:
        nearest_gene_daf = pd.json_normalize(response_data['variantInfo'])
    else:
        nearest_gene_daf = pd.DataFrame({'nearestGene.symbol': [''], 'nearestGeneDistance': [''], 'top_genes': [''], 'mostSevereConsequence': ['']})

    return {'rsid' : rsid,
            'nearest_gene': nearest_gene_daf.squeeze()['nearestGene.symbol'],
            'nearest_gene_distance': nearest_gene_daf.squeeze()['nearestGeneDistance'],
            'top_genes': top_genes_str,
            'most_severe_consequence': nearest_gene_daf.squeeze()['mostSevereConsequence']
            }

# NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints
genes_query = """
query nearest_gene_and_top_genes($variant_id: String!) {
    variantInfo(variantId: $variant_id) {
        nearestGene {
            symbol
        }
        nearestGeneDistance
        mostSevereConsequence
        }
    genesForVariant(variantId: $variant_id) {
        overallScore
        gene {
        symbol
        id
        }
  }
}
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv(snakemake.input[0], sep = '\\s+')

if daf.shape[0] == 0:
    pd.DataFrame(columns = ['rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'nearest_gene', 'nearest_gene_distance', 'top_gene', 'most_severe_consequence']).to_csv(snakemake.output[0], sep = '\t')
else:
    res_dicts = []

    for index, row in daf.iterrows():
         res_dicts.append(query_snp(row.rsid, row.chromosome, row.base_pair_location, row.other_allele, row.effect_allele))

    annot_daf = pd.DataFrame(res_dicts)

    pd.merge(daf, annot_daf, on = 'rsid').to_csv(snakemake.output[0], sep = '\t', index = False)
