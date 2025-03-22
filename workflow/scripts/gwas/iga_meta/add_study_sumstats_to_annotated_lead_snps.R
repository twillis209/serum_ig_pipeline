library(data.table)
setDTthreads(snakemake@threads)
library(stringr)

lead <- fread(snakemake@input$lead)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col

coord_cols <- c(chr_col, bp_col, ref_col, alt_col)

inclusion_wildcards <- snakemake@wildcards[str_detect(names(snakemake@wildcards), '_inclusion')]

studies <- character()

for(x in inclusion_wildcards) {
  if(str_detect(x, 'with_')) {
    studies <- c(studies, (str_split_1(x, '_')[2]))
  }
}

studies <- sort(studies)

cols <- as.vector(outer(c(beta_col, se_col, p_col), studies, FUN = paste, sep = '.'))

cols <- c(coord_cols, cols)

if('epic' %in% studies) {
  cols <- c(cols, 'effect_allele_frequency.epic')
}

if('pietzner' %in% studies) {
  cols <- c(cols, 'effect_allele_frequency.pietzner')
}

if('scepanovic' %in% studies) {
  cols <- c(cols, 'effect_allele_frequency.scepanovic')
}

if('dennis' %in% studies) {
  cols <- c(cols, 'AF1.dennis')
}

if('eldjarn' %in% studies) {
  cols <- c(cols, 'ImpMAF.eldjarn')
}

merged <- fread(snakemake@input$merged, select = cols)

lead_with_sumstats <- merge(lead, merged, by = coord_cols, all.x = T)

names(lead_with_sumstats) <- gsub("effect_allele_frequency|AF1|ImpMAF", "maf", names(lead_with_sumstats))

maf_cols <- names(lead_with_sumstats)[names(lead_with_sumstats) %like% 'maf\\.']
other_cols <- setdiff(names(lead_with_sumstats), maf_cols)
setcolorder(lead_with_sumstats, c(other_cols, sort(maf_cols)))

lead_with_sumstats[, names(.SD) := lapply(.SD, function(x) ifelse(x > 0.5, 1-x, x)), .SDcols = patterns("^maf\\.")]

studies_with_isotype <- paste(studies[studies %in% snakemake@config$studies_with_maf_estimates], snakemake@params$isotype, sep = '-')

sample_sizes <- sapply(studies_with_isotype, function(x) snakemake@config$gwas_datasets[[x]][['samples']], USE.NAMES = T)

lead_with_sumstats[,
  maf.meta := apply(.SD, 1, function(maf_row) {
  # Use the provided sample_sizes vector.
  # Create a weight vector where NA MAF's get a weight of zero.
  effective_weights <- ifelse(is.na(maf_row), 0, sample_sizes)

  # Calculate the weighted average.
  if(sum(effective_weights) > 0) {
    sum(maf_row * sample_sizes, na.rm = TRUE) / sum(effective_weights)
  } else {
    NA_real_
  }
}), .SDcols = patterns("^maf\\.")]


fwrite(lead_with_sumstats, file = snakemake@output[[1]], sep = '\t')
