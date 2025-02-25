import requests
import json
import pandas as pd
import re

gene_query = """
query get_gene_coordinates($gene_name: String!){
    search(queryString: $gene_name) {
        genes {
            id
            symbol
            chromosome
            start
            end
            tss
            fwdStrand
        }
    }
}
"""

base_url = "https://api.genetics.opentargets.org/graphql"

gene_daf = pd.read_csv(snakemake.input[0], sep = '\t', header = 0)

def query_gene(sanitised_gene, tangye_gene):
    print(sanitised_gene)
    variables = {"sanitised_gene": sanitised_gene}

    r = requests.post(base_url, json={"query": gene_query, "variables": variables})

    response_daf = pd.DataFrame(json.loads(r.text)['data']['search']['genes'])
    response_daf['sanitised_gene'] = sanitised_gene
    response_daf['tangye_gene'] = tangye_gene

    return response_daf

def query_gene_with_ensembl(ensembl_id):
    base_url = "https://rest.ensembl.org"
    endpoint = f"/lookup/id/{ensembl_id}?expand=1"

    headers = {
        "Content-Type": "application/json"
    }

    response = requests.get(base_url + endpoint, headers=headers)

    if response.status_code == 200:
        gene_data = response.json()
        if "seq_region_name" in gene_data and "start" in gene_data and "end" in gene_data:
            chromosome = gene_data["seq_region_name"]
            start = gene_data["start"]
            end = gene_data["end"]
            strand = gene_data['strand']
            symbol = gene_data['display_name']
            tss = start if strand == 1 else end

            return pd.DataFrame([{'id': ensembl_id, 'symbol': symbol, 'chromosome': chromosome, 'start': start, 'end': end, 'tss': tss, 'fwdStrand': strand == 1}])
        else:
            return None
    else:
        print("Error:", response.status_code)
        return None

d = []

for index, row in gene_daf.iterrows():
    d.append(query_gene(row.sanitised_gene, row.tangye_gene))

daf = pd.concat(d)

# Handle multiple matches case-by-case
out_daf = daf[daf['sanitised_gene'] == daf['symbol']]

matching_gene_names = out_daf.gene_name

# Handle special cases
ada_id = 'ENSG00000196839'
out_daf = pd.concat([out_daf, query_gene(ada_id, 'ADA').assign(sanitised_gene = 'ADA')])

cd20_id = 'ENSG00000156738'
out_daf = pd.concat([out_daf, query_gene(cd20_id, 'CD20').assign(sanitised_gene = 'CD20')])

# TNFRSF6
fas_id = 'ENSG00000026103'
out_daf = pd.concat([out_daf, query_gene(fas_id, 'TNFRSF6').assign(sanitised_gene = 'TNFRSF6')])

# HMOX aka HMOX1
out_daf = pd.concat([out_daf, daf.query('symbol == \'HMOX1\'')])

# PSEN aka PSEN1
out_daf = pd.concat([out_daf, daf.query('symbol == \'PSEN1\'')])

# C2
c2_id = 'ENSG00000166278'
out_daf = pd.concat([out_daf, query_gene(c2_id, 'C2').assign(sanitised_gene = 'C2')])

# CD3Z, CD247
out_daf = pd.concat([out_daf, query_gene('CD247', 'CD247').assign(sanitised_gene = 'CD3Z')])

# SLP76, LCP2
out_daf = pd.concat([out_daf, query_gene('LCP2', 'LCP2').assign(sanitised_gene = 'SLP76')])

# TRAC
trac_id = 'ENSG00000277734'
out_daf = pd.concat([out_daf, query_gene_with_ensembl(trac_id).assign(sanitised_gene = 'TRAC')])

# NBS1
nbs1_id = 'ENSG00000104320'
out_daf = pd.concat([out_daf, query_gene(nbs1_id, 'NBN').assign(sanitised_gene = 'NBS1')])

pole1_id = 'ENSG00000177084'
out_daf = pd.concat([out_daf, query_gene(pole1_id, 'POLE1').assign(sanitised_gene = 'POLE1')])

rmrp_id = 'ENSG00000277027'
out_daf = pd.concat([out_daf, query_gene_with_ensembl(rmrp_id, 'RMRP').assign(sanitised_gene = 'RMRP')])

## NB: mislabelled in the Tangye sheet as 'ERBB21P' instead of 'ERBB2IP'; the 'I' is for 'interacting'
#erbin_id = 'ENSG00000112851'
#out_daf = pd.concat([out_daf, query_gene(erbin_id, 'ERBIN',).assign(sanitised_gene = 'ERBB2IP')])

cd21_id = 'ENSG00000117322'
out_daf = pd.concat([out_daf, query_gene(cd21_id).assign(sanitised_gene = 'CD21')])

# TNFSF6 is FASLG
tnfsf6_id = 'ENSG00000117560'
out_daf = pd.concat([out_daf, query_gene(tnfsf6_id).assign(sanitised_gene = 'TNFSF6')])

g6pt1_id = 'ENSG00000137700'
out_daf = pd.concat([out_daf, query_gene(g6pt1_id).assign(sanitised_gene = 'G6PT1')])

taz_id = 'ENSG00000018408'
out_daf = pd.concat([out_daf, query_gene(taz_id).assign(sanitised_gene = 'TAZ')])

mkl1_id = 'ENSG00000196588'
out_daf = pd.concat([out_daf, query_gene(mkl1_id).assign(sanitised_gene = 'MKL1')])

# NB: has a trailing space in the name on the Tangye gene list
il12b_id = 'ENSG00000113302'
out_daf = pd.concat([out_daf, query_gene(il12b_id).assign(sanitised_gene = 'IL12B')])

snora31_id = 'ENSG00000199477'
out_daf = pd.concat([out_daf, query_gene_with_ensembl(snora31_id).assign(sanitised_gene = 'SNORA31')])

adar1_id = 'ENSG00000160710'
out_daf = pd.concat([out_daf, query_gene(adar1_id).assign(sanitised_gene = 'ADAR1')])

rnu7_1_id = 'ENSG00000238923'
out_daf = pd.concat([out_daf, query_gene_with_ensembl(rnu7_1_id).assign(sanitised_gene = 'RNU7-1')])

rnu4atac_id = 'ENSG00000264229'
out_daf = pd.concat([out_daf, query_gene_with_ensembl(rnu4atac_id).assign(sanitised_gene = 'RNU4ATAC')])

tmem173_id = 'ENSG00000184584'
out_daf = pd.concat([out_daf, query_gene(tmem173_id).assign(sanitised_gene = 'TMEM173')])

# NB: mislabelled as NICKAPIL instead of NCKAP1L
nckap1l_id = 'ENSG00000123338'
out_daf = pd.concat([out_daf, query_gene(nckap1l_id).assign(sanitised_gene = 'NCKAP1L')])

nola2_id = 'ENSG00000145912'
out_daf = pd.concat([out_daf, query_gene(nola2_id).assign(sanitised_gene = 'NOLA2')])

nola3_id = 'ENSG00000182117'
out_daf = pd.concat([out_daf, query_gene(nola3_id).assign(sanitised_gene = 'NOLA3')])

xrcc9_id = 'ENSG00000221829'
out_daf = pd.concat([out_daf, query_gene(xrcc9_id).assign(sanitised_gene = 'XRCC9')])

terc_id = 'ENSG00000270141'
out_daf = pd.concat([out_daf, query_gene_with_ensembl(terc_id).assign(sanitised_gene = 'TERC')])

out_daf = out_daf.astype({'start': 'int32', 'end': 'int32', 'tss': 'int32'})
out_daf.to_csv(snakemake.output[0], sep = '\t', index = False)
