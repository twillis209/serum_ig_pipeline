library(data.table)
setDTthreads(snakemake@threads)
library(tomics)
library(stringr)

dat <- fread(snakemake@input[[1]])

add_study <- function(study_name) {
  update_meta_analysis_estimates(dat, study_name,
                                 effect_label = snakemake@config$beta_col,
                                 se_label = snakemake@config$se_col,
                                 p_label = snakemake@config$p_col)
}

inclusion_wildcards <- snakemake@wildcards[str_detect(names(snakemake@wildcards), '_inclusion')]

studies <- c()

for(x in inclusion_wildcards) {
  if(str_detect(x, 'with_')) {
    study <- str_split_1(x, '_')[2]
    studies <- c(studies, study)
    add_study(study)
  }
}

studies <- sort(studies)

beta_cols <- paste('beta', studies, sep = '.')
isntna_cols <- paste('isntna', studies, sep = '.')

dat[, (isntna_cols) := lapply(.SD, function(x) !is.na(x)), .SDcols = beta_cols]

studies_with_isotype <- paste(studies, snakemake@params$isotype, sep = '-')

sample_sizes <- sapply(studies_with_isotype, function(x) snakemake@config$gwas_datasets[[x]][['samples']], USE.NAMES = T)

dat[, sample_size := as.matrix(.SD) %*% sample_sizes, .SDcols = isntna_cols]

cols <- c(snakemake@config$chr_col,
          snakemake@config$bp_col,
          snakemake@config$ref_col,
          snakemake@config$alt_col,
          snakemake@config$rsid_col,
          snakemake@config$beta_col,
          snakemake@config$se_col,
          snakemake@config$p_col,
          'sample_size'
          )

fwrite(dat[, ..cols], file = snakemake@output[[1]], sep = '\t')
