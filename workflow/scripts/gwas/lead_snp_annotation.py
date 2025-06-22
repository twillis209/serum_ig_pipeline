import requests
import json
import pandas as pd
import re

"""
I will add:
- genes for which the lead variant is missense;
- genes for which the lead variants is in a qtl credible set;
- if none of the above,  nearest gene
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

    daf = pd.json_normalize(response_data["search"]["hits"]).sort_values(
        by="score", ascending=False
    )

    return daf.loc[0, "id"]

def query_snp_credible_sets(variant_id, size=1000, index=0, include_trans_qtls=False):
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

    r = requests.post(
        base_url, json={"query": credible_sets_query, "variables": variables}
    )

    response_data = json.loads(r.text)["data"]

    if not response_data["variant"]:
        return None
    else:
        daf = pd.json_normalize(response_data["variant"]["qtlCredibleSets"]["rows"])

        daf["locus.rows_extracted"] = daf["locus.rows"].apply(
            lambda x: x[0] if isinstance(x, list) and x else {}
        )

        expanded = daf["locus.rows_extracted"].apply(pd.Series)

        daf = daf.drop(["locus.rows", "locus.rows_extracted"], axis=1).join(expanded)

        daf = daf[daf["isTransQtl"] == include_trans_qtls]

        return daf

def query_snp_missense_and_nearest(variant_id):
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

