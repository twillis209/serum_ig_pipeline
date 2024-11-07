library(data.table)
library(dplyr)
library(httr)

existing_loci <- fread(snakemake@input[['existing']])
new_loci <- fread(snakemake@input[['new']])

bp_flank <- snakemake@params[['window']]/2
novel <- snakemake@params$novel

new_loci[, `:=` (BP.left = BP - bp_flank, BP.right = BP + bp_flank, CHR = as.character(CHR))]
setkey(new_loci, CHR, BP.left, BP.right)
existing_loci[, `:=` (BP.left = BP38 - bp_flank, BP.right = BP38 + bp_flank, CHR = as.character(CHR38))]
setkey(existing_loci, CHR, BP.left, BP.right)

overlaps <- foverlaps(existing_loci, new_loci) %>%
  .[!(CHR == '6' & BP %between% c(24e6, 45e6))] %>%
  .[, .(rsID = i.rsID, topGene, genes, CHR, BP, CHR38, BP38, REF38, ALT38)] %>%
  unique(by = 'rsID') %>%
  .[is.na(topGene), topGene := genes] %>%
  unique(by = 'topGene') %>%
  .[order(CHR38, BP38)]

actually_new_loci <- new_loci[topGene %like% paste(novel, collapse = '|')]

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

for(i in 1:nrow(overlaps)) {
  variables <- list("variantId" = overlaps[i, rsID])
  post_body <- list(query = rsID_coordinate_query, variables = variables)
  r <- POST(url = base_url, body = post_body, encode = 'json')
  res <- content(r)$data$search$variants[[1]]
  overlaps[i, `:=` (CHR19 = res$chromosomeB37, BP19 = res$positionB37)]
}

all_loci <- rbindlist(list(
  overlaps[, .(rsID, topGene, genes, CHR38, BP38, CHR19, BP19, REF38, ALT38)],
  actually_new_loci[, .(rsID, CHR38 = CHR, BP38 = BP, REF38 = REF, ALT38 = ALT, genes = paste(topGene, nearestGene, sep = ', '))]
), fill = T)[order(as.integer(CHR38), BP38)]

fwrite(all_loci, file = snakemake@output[[1]], sep = '\t')
