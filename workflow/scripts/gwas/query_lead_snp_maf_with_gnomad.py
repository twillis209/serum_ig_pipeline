import asyncio
import pandas as pd
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport

daf = pd.read_csv(snakemake.input[0], sep = '\t')

# GraphQL API endpoint
API_URL = "https://gnomad.broadinstitute.org/api"

# Create an HTTP transport
transport = AIOHTTPTransport(url=API_URL)

# Create a GraphQL client
client = Client(transport=transport, fetch_schema_from_transport=True)

# Semaphore to limit concurrent requests
semaphore = asyncio.Semaphore(10)  # 10 requests per minute

async def query_variant_for_nfe_af(rsid):
    query = gql(
        """
        query getVariant($rsId: String!) {
            variant(rsid: $rsId, dataset: gnomad_r4) {
                variant_id
                    reference_genome 
                    genome {
                    af
                    ac
                    an 
                    populations {
                        id
                        ac
                        an
                        homozygote_count
                        hemizygote_count
                    }
                    }
                    }
                }
        """
)

    async with semaphore:  # Ensure only 10 requests per minute
        await asyncio.sleep(6)  # Space requests to stay within the limit (60 sec / 10 reqs = 6 sec per request)
        async with client as session:
            result = await session.execute(query, variable_values={"rsId": rsid})
            nfe_dict = [x for x in result['variant']['genome']['populations'] if x['id'] in ['nfe']][0]
# [{'id': 'nfe', 'ac': 0, 'an': 68016, 'homozygote_count': 0, 'hemizygote_count': 0}]
            af = nfe_dict['ac']/nfe_dict['an']

    return af

    return result

async def main():
    rsids = daf[daf['gnomadNFE'].isna()].rsid.tolist()
    tasks = [query_variant_for_nfe_af(rsid) for rsid in rsids]
    await asyncio.gather(*tasks)
    daf.loc[daf['gnomadNFE'].isna(), 'gnomadNFE'] = tasks

asyncio.run(main())

daf.to_csv(snakemake.output[0], sep = '\t', index = False)
