library(data.table)
library(LDlinkR)

ig <- fread(snakemake@input$ig)

coloc <- fread(snakemake@input$coloc, sep = '\t', header = T)

if(ig[, .N] == 0 | coloc[, .N] == 0) {
  fwrite(data.table(nsnps = integer(), PP.H0.abf = numeric(), PP.H1.abf = numeric(), PP.H2.abf = numeric(), PP.H3.abf = numeric(), PP.H4.abf = numeric(), first_trait = character(), second_trait = character(), ig_snp = character(), non_ig_snp = character(), chromosome = character(), ig_snp_pos = integer(), non_ig_snp_pos = integer(), min_p.first = numeric(), min_p.second = numeric(), ig_snp_effect_ratio = numeric(), non_ig_snp_effect_ratio = numeric(), genes = character()), sep = '\t', file = snakemake@output[[1]])
} else {
  setkey(ig, rsid)

  setkey(coloc, ig_snp)
  coloc[ig, `:=` (genes = i.genes, chromosome = i.chromosome, ig_snp_pos = base_pair_location)]

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
