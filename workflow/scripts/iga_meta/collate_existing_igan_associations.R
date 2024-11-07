library(data.table)
library(httr)

ebi <- fread(snakemake@input[['ebi']])

ebi <- ebi[ , .(region.ebi = REGION, chr.ebi = CHR_ID, bp.ebi = CHR_POS, reported_genes.ebi = `REPORTED GENE(S)`, mapped_genes.ebi = MAPPED_GENE, rsID.ebi = SNPS, risk_allele.ebi = `STRONGEST SNP-RISK ALLELE`, p.ebi = `P-VALUE`)]

ebi <- ebi[p.ebi <= 5e-8]

ebi[rsID.ebi %like% ';', rsID.ebi := tstrsplit(rsID.ebi, split = ';', keep = 1)]

# NB: Based on the sample script provided by Open Targets Genetics here: https://genetics-docs.opentargets.org/data-access/graphql-api#available-endpoints
rsID_coordinate_query = "
query coordinateQuery($variantId: String!) {
    search(queryString: $variantId) {
        variants {
           chromosome
            position
            refAllele
            altAllele
            chromosomeB37
            positionB37
        }
    }
}
"

base_url <- "https://api.genetics.opentargets.org/graphql"

for(i in 1:nrow(ebi)) {
  variables <- list("variantId" = ebi[i, rsID.ebi])
  print(variables)
  post_body <- list(query = rsID_coordinate_query, variables = variables)
  r <- POST(url = base_url, body = post_body, encode = 'json')
  res <- tryCatch(content(r)$data$search$variants[[1]], error = function(e) return(NULL))

  if(!is.null(res)) {
    print(i)
    ebi[i, `:=` (CHR38 = res$chromosome, BP38 = res$position, REF38 = res$refAllele, ALT38 = res$altAllele)]
  }
}

fwrite(ebi[ .(rsID = rsID.ebi, CHR19 = chr.ebi, BP19 = bp.ebi, risk_allele = risk_allele.ebi, CHR38, BP38, REF38, ALT38, reported_genes = reported_genes.ebi, mapped_genes = mapped_genes.ebi, P = p.ebi)], file = snakemake@output[[1]], sep = '\t')
