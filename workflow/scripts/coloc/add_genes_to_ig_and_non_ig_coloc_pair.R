library(data.table)

ig <- fread(snakemake@input$ig)
ig[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(ig)]

non_ig <- fread(snakemake@input$non_ig)
non_ig[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(non_ig)]

coloc <- fread(snakemake@input$coloc)

rbound_genes <- rbindlist(list(ig[, .(rsid, genes)],
                               non_ig[, .(rsid, genes)]))

rbound_genes <- rbound_genes[, .(genes = paste(unique(unlist(strsplit(genes, ", "))), collapse = ", ")), by = rsid][genes != '']

setkey(rbound_genes, rsid)

setkey(coloc, first_snp)
coloc[rbound_genes, genes.first_snp := genes]

setkey(coloc, second_snp)
coloc[rbound_genes, genes.second_snp := genes]

coloc[, max_post := names(.SD)[max.col(.SD, ties.method = 'first')], .SDcols = patterns('PP.H')]

fwrite(coloc, file = snakemake@output[[1]], sep = '\t')
