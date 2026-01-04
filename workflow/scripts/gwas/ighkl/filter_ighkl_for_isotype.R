library(data.table)
setDTthreads(snakemake@threads)

isotype <- snakemake@params[['isotype']]

dat <- fread(snakemake@input[[1]], sep = '\t')

base_cols <- c('rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele')

suffix <- paste0('.', isotype, '_meta')

cols_to_keep <- c(
  base_cols,
  paste0('beta', suffix),
  paste0('standard_error', suffix),
  paste0('p_value', suffix)
)

dat <- dat[, ..cols_to_keep]

setnames(dat, 
         old = c(paste0('beta', suffix), paste0('standard_error', suffix), paste0('p_value', suffix)),
         new = c('beta', 'standard_error', 'p_value'))

dat <- dat[!is.na(p_value)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
