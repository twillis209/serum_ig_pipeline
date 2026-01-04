library(data.table)
setDTthreads(snakemake@threads)
library(stringr)

dat <- fread(snakemake@input[[1]], sep = '\t')

base_cols <- c('rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele')

study <- snakemake@wildcards$study

# Non-meta-analysis dataset
if (!str_detect(study, "meta")) {
  suffix <- gsub("_", ".", study)
  study_with_dash <- gsub("_", "-", study)
  dat[, sample_size := snakemake@config$gwas_datasets[[study_with_dash]][["samples"]]]
  setnames(dat, "sample_size", paste0("sample_size", suffix))
}

cols_to_keep <- c(
  base_cols,
  paste0('beta', suffix),
  paste0('standard_error', suffix),
  paste0("p_value", suffix),
  paste0("sample_size", suffix)
)

dat <- dat[, ..cols_to_keep]

setnames(dat, 
         old = c(paste0('beta', suffix), paste0('standard_error', suffix), paste0('p_value', suffix), paste0('sample_size', suffix)),
         new = c('beta', 'standard_error', 'p_value', 'sample_size')
         )

dat <- dat[!is.na(p_value)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
