library(data.table)
library(LDlinkR)

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

for(i in seq_len(coloc[, .N])) {
  r2 <- tryCatch({
  ld <- LDpop(
    var1 = coloc[i, first_snp],
    var2 = coloc[i, second_snp],
    pop = "EUR",
    r2d = "r2",
    token = snakemake@config$LDlink$token,
    file = FALSE,
    genome_build = "grch38_high_coverage"
  )
  ld[ld$Abbrev == 'EUR', 'R2']
  },
  error = function(e) { message(e); NA })

 set(coloc, i, "r2", r2)
}

fwrite(coloc, file = snakemake@output[[1]], sep = '\t')
