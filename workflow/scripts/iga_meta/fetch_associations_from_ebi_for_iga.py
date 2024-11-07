import requests
import json
import pandas as pd

# TODO findByEfoTrait

base_url = "https://www.ebi.ac.uk/gwas/rest/api/associations/search/findByEfoTrait"
efoTrait = "EFO_0004912"

# Send an HTTP GET request to retrieve study associations
response = requests.get(f"{base_url}?{efoTrait}")

# Check if the request was successful (HTTP status code 200)
if response.status_code == 200:
    # Parse the JSON response
    data = response.json()

    # Access the associations data in the response
    associations = data['_embedded']['associations']

    with open(snakemake.output.json, 'w') as fh:
        fh.write(json.dumps(associations))

    entries = []

    # Process and print the associations
    for association in associations:
        # Access association data fields as needed
        gene_names = []

        for gene in association['loci'][0]['authorReportedGenes']:
            gene_names.append(gene['geneName'])

        gene_names = ','.join(gene_names)

        entries.append({
            'riskFrequency': association['riskFrequency'],
            'pvalue': association['pvalueMantissa']*pow(10, association['pvalueExponent']),
            'standardError': association['standardError'],
            'range': association['range'],
            'orPerCopyNum': association['orPerCopyNum'],
            'betaNum': association['betaNum'],
            'betaUnit': association['betaUnit'],
            'betaDirection': association['betaDirection'],
            'rsID' : association['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[0],
            'risk_allele' : association['loci'][0]['strongestRiskAlleles'][0]['riskAlleleName'].split('-')[1],
            'risk_frequency' : association['loci'][0]['strongestRiskAlleles'][0]['riskFrequency'],
            'gene_names' : gene_names
        })
    daf = pd.DataFrame(entries)
    daf.to_csv(snakemake.output.tsv, sep = '\t', index = False)
else:
    print(f"Error: Unable to retrieve data. Status code: {response.status_code}")

# Close the response
response.close()
