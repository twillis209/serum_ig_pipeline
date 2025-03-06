library(data.table)

all_pairs <- fread(snakemake@input$all_pairs)
setkey(all_pairs, first_snp)

iga <- fread(snakemake@input$iga)
iga[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(iga)]
setkey(iga, rsid)
igg <- fread(snakemake@input$igg)
igg[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(igg)]
setkey(igg, rsid)
igm <- fread(snakemake@input$igm)
igm[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(igm)]
setkey(igm, rsid)

