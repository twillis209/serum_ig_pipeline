import requests
import json
import pandas as pd
import re

chr_col = snakemake.config['chr_col']
bp_col = snakemake.config['bp_col']
ref_col = snakemake.config['ref_col']
alt_col = snakemake.config['alt_col']
beta_col = snakemake.config['beta_col']
se_col = snakemake.config['se_col']
p_col = snakemake.config['p_col']
rsid_col = snakemake.config['rsid_col']
snp_col = 'variant_id'

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

daf = pd.read_csv(snakemake.input[0], sep = '\\s+')

if daf.shape[0] == 0:
    pd.DataFrame(columns = ['variant_id', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'rsid', 'nearestGene', 'nearestGeneDistance', 'topGene', 'topGeneId', 'mostSevereConsequence', 'gnomadFIN', 'gnomadNFE', 'gnomadAFR', 'gnomadAMR', 'gnomadASJ', 'gnomadEAS', 'gnomadNFEEST', 'gnomadNFENWE', 'gnomadNFESEU']).to_csv(snakemake.output[0], sep = '\t')
else:
    daf[snp_col] = daf[[chr_col, bp_col, ref_col, alt_col]].astype(str).agg('_'.join, axis=1)

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
        res_dict = query_snp(row[snp_col], row[chr_col], row[bp_col], row[ref_col], row[alt_col])

        if res_dict['variantInfo'] is None:
            res_dict = query_snp(row.SNP, row[chr_col], row[bp_col], row[ref_col], row[alt_col], ref_alt_order = False)

        result_dict[row[snp_col]] = res_dict

#        if rsid_col in row.keys():
#            result_dict[row[snp_col]]['variantInfo'][rsid_col] = row[rsid_col]

        result_dict[row[snp_col]]['variantInfo'][chr_col] = row[chr_col]
        result_dict[row[snp_col]]['variantInfo'][bp_col] = row[bp_col]
        result_dict[row[snp_col]]['variantInfo'][ref_col] = row[ref_col]
        result_dict[row[snp_col]]['variantInfo'][alt_col] = row[alt_col]
        result_dict[row[snp_col]]['variantInfo'][rsid_col] = row[rsid_col]
        result_dict[row[snp_col]]['variantInfo'][beta_col] = row[beta_col]
        result_dict[row[snp_col]]['variantInfo'][se_col] = row[se_col]
        result_dict[row[snp_col]]['variantInfo'][p_col] = row[p_col]

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


        d.append(
            {
                snp_col: k,
                rsid_col: v['variantInfo'][rsid_col],
                chr_col: v['variantInfo'][chr_col],
                bp_col: v['variantInfo'][bp_col],
                ref_col: v['variantInfo'][ref_col],
                alt_col: v['variantInfo'][alt_col],
                beta_col: v['variantInfo'][beta_col],
                se_col: v['variantInfo'][se_col],
                p_col: v['variantInfo'][p_col],
                'nearestGene':  nearest_gene,
                'nearestGeneDistance' : v['variantInfo']['nearestGeneDistance'],
                'topGene' : top_gene,
                'topGeneId' : top_gene_id,
                'mostSevereConsequence' : v['variantInfo']['mostSevereConsequence']
            }
            )

    pd.DataFrame(d).to_csv(snakemake.output[0], sep = '\t', index = False)
