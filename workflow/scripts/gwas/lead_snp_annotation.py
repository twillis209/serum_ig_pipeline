import requests
import json
import pandas as pd
import re

chr_col = snakemake.config.chr_col
bp_col = snakemake.config.bp_col
snp_col = snakemake.config.id_col
ref_col = snakemake.config.ref_col
alt_col = snakemake.config.alt_col

# NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints

index_variants_and_studies_query = """
    query annotateLeadSnp($inputVariantId: String!) {
        indexVariantsAndStudiesForTagVariant(variantId: $inputVariantId) {
            associations {
            indexVariant {
            nearestGene {
                id
                symbol
            }
            mostSevereConsequence
            id
            }
            study {
                source
                traitReported
                pmid
                pubDate
                pubTitle
                pubAuthor
                hasSumstats
                nInitial
                nReplication
                nCases
                traitCategory
                numAssocLoci
            }
            pval
            beta
        }
  }
}
"""

variant_query = """
    query annotateLeadSnp($inputVariantId: String!) {
        variantInfo(variantId: $inputVariantId) {
            rsId
            nearestGene {
                id
                symbol
                bioType
            }
            nearestGeneDistance
            mostSevereConsequence
            gnomadAFR
            gnomadAMR
            gnomadASJ
            gnomadEAS
            gnomadFIN
            gnomadNFE
            gnomadNFEEST
            gnomadNFENWE
            gnomadNFESEU
            }
}
"""

genes_for_variant_query = """
query annotateLeadSnp($inputVariantId: String!){
    genesForVariant(variantId: $inputVariantId) {
        overallScore
        gene {
        id
        symbol
        bioType
        description
        }
        functionalPredictions {
        typeId
        sourceId
        aggregatedScore
        tissues {
            tissue {
            id
            }
            maxEffectLabel
            maxEffectScore
        }
        }
        distances {
        sourceId
        aggregatedScore
        tissues {
            tissue {
            id
            name
            }        
            distance 
        }
        }
        intervals {
        typeId
        sourceId
        tissues {
            score
        }
        }
        qtls {
        sourceId
        typeId
        aggregatedScore
                tissues {
            tissue {
            id
            }
                beta
                pval
                }
        }
    }
    }
"""

base_url = "https://api.genetics.opentargets.org/graphql"

daf = pd.read_csv(snakemake.input[0], delim_whitespace = True)

daf.rename(columns = {chr_col: 'CHR', bp_col: 'BP', snp_col: 'SNP', ref_col: 'REF', alt_col: 'ALT'}, inplace = True)

result_dict = {}

def query_snp(snp_id, chrom, bp, ref, alt, ref_alt_order = True):
    if ref_alt_order:
        variables = {"inputVariantId": '_'.join([str(x) for x in [chrom, bp, ref, alt]])}
    else:
        variables = {"inputVariantId": '_'.join([str(x) for x in [chrom, bp, alt, ref]])}

    r = requests.post(base_url, json={"query": variant_query, "variables": variables})

    variant_response_data = json.loads(r.text)['data']['variantInfo']

    r = requests.post(base_url, json={"query": index_variants_and_studies_query, "variables": variables})

    index_variants_and_studies_response_data = json.loads(r.text)['data']['indexVariantsAndStudiesForTagVariant']

    r = requests.post(base_url, json={"query": genes_for_variant_query, "variables": variables})

    genes_for_variant_response_data = json.loads(r.text)['data']['genesForVariant']

    return {'variantInfo' : variant_response_data, 'indexVariantsAndStudiesForTagVariant': index_variants_and_studies_response_data, 'genesForVariant': genes_for_variant_response_data}

for index, row in daf.iterrows():
    res_dict = query_snp(row.SNP, row.CHR, row.BP, row.REF, row.ALT)

    if res_dict['variantInfo'] is None or res_dict['variantInfo']['rsId'] == 'null':
        res_dict = query_snp(row.SNP, row.CHR, row.BP, row.REF, row.ALT, ref_alt_order = False)

    result_dict[row.SNP] = res_dict

with open(snakemake.output.annotations, 'w') as f:
    json.dump(result_dict, f)

d = []

for k,v in result_dict.items():
    study_string = ';'.join([f"{x['study']['traitReported']}" for x in v['indexVariantsAndStudiesForTagVariant']['associations'] if x['pval'] <= 5e-8])
    study_string = study_string.replace('"', '')
    study_string = f"\"{study_string}\""

    nearest_gene = '' if v['variantInfo']['nearestGene'] is None else v['variantInfo']['nearestGene']['symbol']

    top_gene = ""
    top_gene_id = ""
    top_score = 0

    if v['genesForVariant']:
        for x in v['genesForVariant']:
            if x['overallScore'] > top_score:
                top_score = x['overallScore']
                top_gene = x['gene']['symbol']
                top_gene_id = x['gene']['id']


    # TODO handle case where v['variantInfo'] is None
    d.append(
        {
            'SNPID': k,
            'rsID': v['variantInfo']['rsId'],
            'nearestGene':  nearest_gene,
            'nearestGeneDistance' : v['variantInfo']['nearestGeneDistance'],
            'topGene' : top_gene,
            'topGeneId' : top_gene_id,
            'mostSevereConsequence' : v['variantInfo']['mostSevereConsequence'],
            'gnomadFIN' : v['variantInfo']['gnomadFIN'],
            'gnomadNFE' : v['variantInfo']['gnomadNFE'],
            'gnomadAFR' : v['variantInfo']['gnomadAFR'],
            'gnomadAMR' : v['variantInfo']['gnomadAMR'],
            'gnomadASJ' : v['variantInfo']['gnomadASJ'],
            'gnomadEAS' : v['variantInfo']['gnomadEAS'],
            'gnomadNFEEST' : v['variantInfo']['gnomadNFEEST'],
            'gnomadNFENWE' : v['variantInfo']['gnomadNFENWE'],
            'gnomadNFESEU' : v['variantInfo']['gnomadNFESEU'],
            'associatedTraits' : study_string
        }
        )

meta_daf = pd.DataFrame(d)

meta_daf = meta_daf.merge(right = daf, how = 'right', left_on = 'SNPID', right_on = 'SNP')

meta_daf.to_csv(snakemake.output.rsIDs, sep = '\t', index = False)
