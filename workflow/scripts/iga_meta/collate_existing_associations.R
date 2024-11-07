library(data.table)
library(httr)

ebi <- fread(snakemake@input[['ebi']])
liu <- fread(snakemake@input[['liu']])
pietzner <- fread(snakemake@input[['pietzner']])

flank <- 1e5

ebi <- rbindlist(list(ebi, pietzner))

ebi <- ebi[, .(region.ebi = REGION, chr.ebi = CHR_ID, bp.ebi = CHR_POS, reported_genes.ebi = `REPORTED GENE(S)`, mapped_genes.ebi = MAPPED_GENE, rsID.ebi = SNPS, risk_allele.ebi = `STRONGEST SNP-RISK ALLELE`, beta.ebi = `OR or BETA`, p.ebi = `P-VALUE`)]
ebi[, risk_allele.ebi := tstrsplit(risk_allele.ebi, split = '-', keep = 2)]
ebi[, `:=` (bp.left.ebi = bp.ebi - flank, bp.right.ebi = bp.ebi + flank)]
ebi <- ebi[p.ebi <= 5e-8]
setkey(ebi, chr.ebi, bp.left.ebi, bp.right.ebi)

liu <- liu[, .(chr.liu = as.character(CHR), bp.liu = `BP (hg19)`, rsID.liu = SNP, effect_allele.liu = `Effect Allele`, p.liu = `P-value FE`, reported_genes.liu = Locus, beta.liu = Beta)]
liu[, `:=` (bp.left.liu = bp.liu - flank, bp.right.liu = bp.liu + flank)]
setkey(liu, chr.liu, bp.left.liu, bp.right.liu)

overlaps <- foverlaps(ebi, liu, mult = 'first')

liu_matches <- overlaps[!is.na(rsID.liu)]

join <- rbindlist(list(overlaps, liu[!(rsID.liu %in% liu_matches$rsID.liu)]), fill = T)
join[, chr := ifelse(is.na(chr.liu), chr.ebi, chr.liu)]
join[, bp := ifelse(is.na(bp.liu), bp.ebi, bp.liu)]
join[, chr := factor(chr, levels = c(1:22, 'X'))]
# For the two rows where we have information from both, just take Liu
join[!is.na(bp.liu) & !is.na(bp.ebi), `:=` (chr19 = chr.liu, bp19 = bp.liu, rsID = rsID.liu, genes = reported_genes.liu)]
join[is.na(chr19) & is.na(bp19) & is.na(bp.ebi), `:=` (chr19 = chr.liu, bp19 = bp.liu)]
join[is.na(chr19) & is.na(bp19) & is.na(bp.liu), `:=` (chr19 = chr.ebi, bp19 = bp.ebi)]
join[is.na(genes) & !is.na(reported_genes.liu), genes := reported_genes.liu]
join[is.na(genes) & !is.na(mapped_genes.ebi), genes := mapped_genes.ebi]
join[is.na(rsID) & !is.na(rsID.ebi), rsID := rsID.ebi]
join[is.na(rsID) & !is.na(rsID.liu), rsID := rsID.liu]
join[, c("bp.left.liu", "bp.right.liu", 'reported_genes.ebi', 'bp.left.ebi', 'bp.right.ebi', 'chr.ebi', 'chr.liu') := NULL]
# Drop repeat entries where there are multiple reported risk/effect alleles
join <- unique(join[, .(rsID, chr19, bp19, effect_allele.liu, risk_allele.ebi, genes, beta.ebi, beta.liu, p.ebi, p.liu)][order(chr19, bp19)], by = c('chr19', 'bp19'))

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

for(i in 1:nrow(join)) {
  variables <- list("variantId" = join[i, rsID])
  post_body <- list(query = rsID_coordinate_query, variables = variables)
  r <- POST(url = base_url, body = post_body, encode = 'json')
  res <- content(r)$data$search$variants[[1]]
  join[i, `:=` (chr38 = res$chromosome, bp38 = res$position, ref38 = res$refAllele, alt38 = res$altAllele, chr19.otg = res$chromosomeB37, bp19.otg = res$positionB37)]
}

join[, chr38 := factor(chr38, levels = c(1:22, 'X'))]

fwrite(join[order(chr38, bp38), .(rsID, CHR19 = chr19.otg, BP19 = bp19.otg, CHR38 = chr38, BP38 = bp38, REF38 = ref38, ALT38 = alt38, risk_allele.ebi, effect_allele.liu, genes, p.ebi, p.liu, beta.ebi, beta.liu)], file = snakemake@output[[1]], sep = '\t')
