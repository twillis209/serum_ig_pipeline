library(data.table)
setDTthreads(snakemake@threads)
library(stringr)

inclusion_wildcards <- snakemake@wildcards[str_detect(names(snakemake@wildcards), '_inclusion')]

studies <- c()

coord_cols <- c(snakemake@config$chr_col,
          snakemake@config$bp_col,
          snakemake@config$ref_col,
          snakemake@config$alt_col,
          snakemake@config$rsid_col
          )

for(x in inclusion_wildcards) {
  if(str_detect(x, 'with_')) {
    studies <- c(studies, str_split_1(x, '_')[2])
  }
}

studies <- sort(studies)

beta_cols <- paste('beta', studies, sep = '.')

cols <- c(coord_cols, beta_cols)

dat <- fread(snakemake@input[[1]], select = cols)

dat[, (beta_cols) := lapply(.SD, function(x) !is.na(x)), .SDcols = beta_cols]

studies_with_isotype <- paste(studies, snakemake@params$isotype, sep = '-')

sample_sizes <- sapply(studies_with_isotype, function(x) snakemake@config$gwas_datasets[[x]][['samples']], USE.NAMES = T)

dat[, 'sample_size' := as.matrix(.SD) %*% sample_sizes, .SDcols = beta_cols]

out_cols <- c(coord_cols, 'sample_size')

fwrite(dat[, ..out_cols], sep = '\t', file = snakemake@output[[1]])

