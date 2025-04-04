library(data.table)
library(LDlinkR)

ig <- fread(snakemake@input$ig)
ig[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(ig)]

coloc <- fread(snakemake@input$coloc)

setkey(ig, rsid)

setkey(coloc, ig_snp)
coloc[ig, genes := genes]

coloc[, max_post := names(.SD)[max.col(.SD, ties.method = 'first')], .SDcols = patterns('PP.H')]

for(i in seq_len(coloc[, .N])) {
  r2 <- tryCatch({
  ld <- LDpop(
    var1 = coloc[i, ig_snp],
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
