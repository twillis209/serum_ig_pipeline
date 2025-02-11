import asyncio
import pandas as pd
from numpy import nan
from gql import Client, gql
from gql.transport.aiohttp import AIOHTTPTransport
from gql.transport.exceptions import TransportQueryError

daf = pd.read_csv(snakemake.input[0], sep = '\t')

# GraphQL API endpoint
API_URL = "https://gnomad.broadinstitute.org/api"

# Create an HTTP transport
transport = AIOHTTPTransport(url=API_URL)

# Create a GraphQL client
client = Client(transport=transport, fetch_schema_from_transport=True)

# Semaphore to limit concurrent requests
semaphore = asyncio.Semaphore(10)  # 10 requests per minute

query_string = """
        query getVariant($rsId: String!) {
            variant(rsid: $rsId, dataset: gnomad_r4) {
                variant_id
                rsid
                chrom
                pos
                ref
                alt
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

async def query_variant_for_nfe_af(session, rsid):
    query = gql(query_string)
    async with semaphore:
        await asyncio.sleep(6)  # Space requests to stay within the limit (60 sec / 10 reqs = 6 sec per request)
        try:
            result = await session.execute(query, variable_values={"rsId": rsid})
            nfe_dict = [x for x in result['variant']['genome']['populations'] if x['id'] in ['nfe']][0]
            af = nfe_dict['ac']/nfe_dict['an']
            return pd.DataFrame([{'variant_id': result['variant']['variant_id'], 'rsid': result['variant']['rsid'], 'chrom': result['variant']['chrom'], 'pos': result['variant']['pos'], 'ref': result['variant']['ref'], 'alt': result['variant']['alt'], 'af': af}])
        except TransportQueryError:
            return nan

async def main():
    async with client as session:
        rsids = daf.rsid.tolist()
        tasks = [query_variant_for_nfe_af(session, rsid) for rsid in rsids]
        results = await asyncio.gather(*tasks)

        return results

dafs = asyncio.run(main())

res_daf = pd.concat(dafs).to_csv(snakemake.output[0], sep = '\t', index = False)
