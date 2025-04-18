library(data.table)
library(LDlinkR)

ig <- fread(snakemake@input$ig)

if(ig[, .N] == 0) {
  fwrite(data.table(nsnps = integer(), PP.H0.abf = numeric(), PP.H1.abf = numeric(), PP.H2.abf = numeric(), PP.H3.abf = numeric(), PP.H4.abf = numeric(), first_trait = character(), second_trait = character(), ig_snp = character(), non_ig_snp = character(), min_p.first = numeric(), min_p.second = numeric()), sep = '\t', file = snakemake@output[[1]])
} else {
  ig[, genes := paste(unique(sort(c(topGene, nearestGene))), collapse = ', '), by = 1:nrow(ig)]

  coloc <- fread(snakemake@input$coloc, sep = '\t', header = T)

  setkey(ig, rsid)

  setkey(coloc, ig_snp)
  coloc[ig, genes := genes]

  coloc[, max_post := names(.SD)[max.col(.SD, ties.method = 'first')], .SDcols = patterns('PP.H')]

  for(i in seq_len(coloc[, .N])) {
    r2 <- tryCatch({
    ld <- LDpop(
      var1 = coloc[i, ig_snp],
      var2 = coloc[i, non_ig_snp],
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
}
