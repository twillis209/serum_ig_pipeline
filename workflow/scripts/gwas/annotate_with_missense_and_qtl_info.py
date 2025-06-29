import requests
import json
import pandas as pd
import re
import asyncio
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport
from gql.transport.exceptions import TransportQueryError

"""
I will add:
- genes for which the lead variant is missense;
- genes for which the lead variants is in a qtl credible set;
- if none of the above,  nearest gene
"""
base_url = "https://api.platform.opentargets.org/api/v4/graphql"


def find_variant_id(row):
    rsid = row["rsid"]

    variant_id = query_snp_id_by_rsid(rsid)

    if variant_id:
        print(f"Variant ID {variant_id} found for rsid: {rsid}")
    # No rsid
    else:
        print(f"Variant ID not found for rsid: {rsid}, concatenating...")
        variant_id = "_".join(
            row[
                ["chromosome", "base_pair_location", "other_allele", "effect_allele"]
            ].astype(str)
        )

        variant_id = query_snp_id_by_variant_id(variant_id)

        # other_effect order doesn't work
        if not variant_id:
            variant_id = "_".join(
                row[
                    [
                        "chromosome",
                        "base_pair_location",
                        "effect_allele",
                        "other_allele",
                    ]
                ].astype(str)
            )

            variant_id = query_snp_id_by_variant_id(variant_id)

            # No luck
            if not variant_id:
                print(f"Variant ID not found for rsid: {rsid}, skipping...")

    return variant_id


def query_snp_id_by_variant_id(variant_id):
    variant_id_query = """
    query variant_id_query($variant_id: String!) {
        search(queryString: $variant_id, entityNames: ["variant"]) {
            hits {
                    id
                    prefixes
                    score
                }
        }
    }
    """

    variables = {"variant_id": variant_id}

    r = requests.post(
        base_url, json={"query": variant_id_query, "variables": variables}
    )

    response_data = json.loads(r.text)["data"]

    # return response_data

    if (
        response_data["search"]["hits"] is None
        or len(response_data["search"]["hits"]) == 0
    ):
        return None

    daf = pd.json_normalize(
        response_data["search"]["hits"], record_path="prefixes", meta=["id", "score"]
    ).sort_values(by="score", ascending=False)

    daf.rename(columns={0: "prefix"}, inplace=True)

    if daf[daf["prefix"] == variant_id].shape[0] == 1:
        return daf[daf["prefix"] == variant_id].iloc[0]["id"]
    elif daf[daf["prefix"] == variant_id].shape[0] > 0:
        # If there are multiple matches, return the one with the highest score
        daf = daf[daf["prefix"] == variant_id].sort_values(by="score", ascending=False)

        return daf.iloc[0]["id"]
    else:
        return None


def query_snp_id_by_rsid(rsid):
    rsid_query = """
    query variant_id_query($rsid: String!) {
        search(queryString: $rsid, entityNames: ["variant"]) {
            hits {
                    id
                    prefixes
                    score
                }
        }
    }
    """

    variables = {"rsid": rsid}

    r = requests.post(base_url, json={"query": rsid_query, "variables": variables})

    response_data = json.loads(r.text)["data"]

    # return response_data

    if (
        response_data["search"]["hits"] is None
        or len(response_data["search"]["hits"]) == 0
    ):
        return None

    daf = pd.json_normalize(
        response_data["search"]["hits"], record_path="prefixes", meta=["id", "score"]
    ).sort_values(by="score", ascending=False)

    daf.rename(columns={0: "prefix"}, inplace=True)

    if daf[daf["prefix"] == rsid].shape[0] == 1:
        return daf[daf["prefix"] == rsid].iloc[0]["id"]
    elif daf[daf["prefix"] == rsid].shape[0] > 0:
        # If there are multiple matches, return the one with the highest score
        daf = daf[daf["prefix"] == rsid].sort_values(by="score", ascending=False)

        return daf.iloc[0]["id"]
    else:
        return None


def query_snp_credible_sets(variant_id, size=1000, index=0, include_trans_qtls=False):
    credible_sets_query = """
    query QTLCredibleSetsQuery($variant_id: String!, $size: Int!, $index: Int!) {
        variant(variantId: $variant_id) {
        id
        qtlCredibleSets: credibleSets(
        studyTypes: [scsqtl, sceqtl, scpqtl, sctuqtl, sqtl, eqtl, pqtl, tuqtl],
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

    r = requests.post(
        base_url, json={"query": credible_sets_query, "variables": variables}
    )

    response_data = json.loads(r.text)["data"]

    if (
        response_data["variant"] is None
        or response_data["variant"]["qtlCredibleSets"]["count"] == 0
    ):
        return None
    else:
        cs_daf = pd.json_normalize(response_data["variant"]["qtlCredibleSets"]["rows"])

        cs_daf["locus.rows_extracted"] = cs_daf["locus.rows"].apply(
            lambda x: x[0] if isinstance(x, list) and x else {}
        )

        expanded = cs_daf["locus.rows_extracted"].apply(pd.Series)

        cs_daf = cs_daf.drop(["locus.rows", "locus.rows_extracted"], axis=1).join(
            expanded
        )

        cs_daf = cs_daf[cs_daf["isTransQtl"] == include_trans_qtls]

        return cs_daf


def query_snp_missense_and_nearest_gene(variant_id):
    missense_and_nearest_gene_query = """
    query missense_and_nearest_gene($variant_id: String!) {
        variant(variantId: $variant_id) {
            id
            mostSevereConsequence {
            label
            id
            }
            transcriptConsequences {
            target {

                approvedSymbol
            }
            variantConsequences {
                # id
                label
            }
            distanceFromFootprint
            distanceFromTss
            target {
                id
                approvedSymbol
                biotype
            }
            }
        }
        }
    """

    variables = {"variant_id": variant_id}

    r = requests.post(
        base_url,
        json={"query": missense_and_nearest_gene_query, "variables": variables},
    )

    response_data = json.loads(r.text)["data"]

    daf = pd.json_normalize(
        response_data["variant"], record_path=["transcriptConsequences"]
    )

    daf["variantConsequences_extracted"] = daf["variantConsequences"].apply(
        lambda x: x[0] if isinstance(x, list) and x else {}
    )

    expanded = daf["variantConsequences_extracted"].apply(pd.Series)

    daf = daf.drop(
        ["variantConsequences", "variantConsequences_extracted"], axis=1
    ).join(expanded)

    daf.rename(
        columns={
            "label": "variantConsequence",
        },
        inplace=True,
    )

    daf.sort_values("distanceFromFootprint", inplace=True)

    return daf


def fetch_genes_using_open_targets_api(daf):
    daf["missense_gene"] = ""
    daf["qtl_genes"] = ""

    reason_dicts = []

    for index, row in daf.iterrows():
        variant_id = find_variant_id(row)

        if variant_id is None:
            print(f"Variant ID not found for row {index}, skipping...")
            continue

        gene = []

        credible_sets = query_snp_credible_sets(variant_id)

        if credible_sets is not None:
            credible_sets = credible_sets[credible_sets["is95CredibleSet"] == True]

        missense_and_nearest = query_snp_missense_and_nearest_gene(variant_id)

        if (
            missense_and_nearest[
                missense_and_nearest["variantConsequence"] == "missense_variant"
            ].shape[0]
            > 0
        ):
            missense_genes = missense_and_nearest[
                missense_and_nearest["variantConsequence"] == "missense_variant"
            ]["target.approvedSymbol"].tolist()


            daf.loc[index, "missense_gene"] = ",".join(list(set(missense_genes)))

            reason_dicts.append(
                {"rsid": row["rsid"], "gene": missense_genes, "reason": "missense"}
            )

        if credible_sets is not None:
            cs_genes = credible_sets["study.target.approvedSymbol"].tolist()

            daf.loc[index, "qtl_genes"] = ",".join(list(set(cs_genes)))

            reason_dicts.append({"rsid": row["rsid"], "gene": cs_genes, "reason": "cs"})

        # if len(gene) == 0:
        #     nearest_gene = missense_and_nearest.head(1)[
        #         "target.approvedSymbol"
        #     ].tolist()

        #     gene += nearest_gene

        #     reason_dicts.append(
        #         {"rsid": row["rsid"], "gene": nearest_gene, "reason": "nearest"}
        #     )

    # reason_daf = pd.DataFrame(reason_dicts)

    # reason_daf = reason_daf.explode("gene")

    return daf

daf = pd.read_csv(snakemake.input[0], sep="\\s+")

gene_daf = fetch_genes_using_open_targets_api(daf)

gene_daf.to_csv(snakemake.output[0], sep = '\t', index = False)

# def query_nearest_genes_using_gnomad_api(rsids):
#     gnomad_url = "https://gnomad.broadinstitute.org/api"
#     transport = AIOHTTPTransport(url=gnomad_url, ssl=False)
#     client = Client(transport=transport, fetch_schema_from_transport=True)
#     semaphore = asyncio.Semaphore(5)

#     async def fetch_variant(rsid, session, semaphore):
#         query_str = """
#             query genes($rsid: String!) {
#               variant(rsid: $rsid, dataset: gnomad_r4) {
#                 variant_id
#                 rsid
#                 chrom
#                 pos
#                 ref
#                 alt
#                 sortedTranscriptConsequences {
#                   consequence_terms
#                   gene_id
#                   gene_symbol
#                 }
#               }
#             }
#         """
#         query = gql(query_str)
#         async with semaphore:
#             try:
#                 result = await session.execute(query, variable_values={"rsid": rsid})
#                 return result if result is not None else None
#             except TransportQueryError:
#                 return None

#     async def main():
#         async with client as session:
#             tasks = [fetch_variant(rsid, session, semaphore) for rsid in rsids]
#             results = await asyncio.gather(*tasks)
#             return results

#     return asyncio.run(main())
