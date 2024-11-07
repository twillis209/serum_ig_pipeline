library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
snp_col <- snakemake@params[['snp_col']]
p_col <- snakemake@params[['p_col']]
beta_col <- snakemake@params[['beta_col']]
se_col <- snakemake@params[['se_col']]

distance_window <- snakemake@params[['distance_window']]
index_threshold <- snakemake@params[['index_threshold']]

dat <- fread(snakemake@input[['gwas']], sep = '\t', select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col, snp_col))

dat[, (chr_col) := as.character(get(chr_col))]
dat[, (bp_col) := as.integer(get(bp_col))]

dat <- dat[!(get(snp_col) %in% snakemake@params$snps_to_ignore)]

dat <- dat[get(p_col) <= index_threshold]

dat <- dat[order(get(p_col))]

for(i in unique(dat[[chr_col]])) {
  j <- 1

  while(j <= dat[get(chr_col) == i, .N]) {
    bp <- dat[get(chr_col) == i][j][[bp_col]]

    dat <- dat[!(get(chr_col) == i & get(bp_col) != bp & get(bp_col) %between% c(bp-distance_window/2, bp+distance_window/2))]

    j <- j + 1
  }
}

dat <- dat[order(get(chr_col), get(bp_col))]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
