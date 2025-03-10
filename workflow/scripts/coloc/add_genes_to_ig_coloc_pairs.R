library(data.table)

all_pairs <- fread(snakemake@input$all_pairs)

iga <- fread(snakemake@input$iga)
iga[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(iga)]
igg <- fread(snakemake@input$igg)
igg[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(igg)]
igm <- fread(snakemake@input$igm)
igm[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(igm)]

rbound_genes <- rbindlist(list(iga[, .(rsid, genes)],
                               igg[, .(rsid, genes)],
                               igm[, .(rsid, genes)]))

rbound_genes <- rbound_genes[, .(genes = paste(unique(unlist(strsplit(genes, ", "))), collapse = ", ")), by = rsid][genes != '']

setkey(rbound_genes, rsid)

setkey(all_pairs, first_snp)
all_pairs[rbound_genes, genes.first_snp := genes]

setkey(all_pairs, second_snp)
all_pairs[rbound_genes, genes.second_snp := genes]

all_pairs[, max_post := names(.SD)[max.col(.SD, ties.method = 'first')], .SDcols = patterns('PP.H')]
all_pairs[, max_post := gsub('PP\\.|\\.abf', '', max_post)]

fwrite(all_pairs, file = snakemake@output[[1]], sep = '\t')
