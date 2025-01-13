library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
p_col <- snakemake@config$p_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col

distance_window <- snakemake@params[['distance_window']]
index_threshold <- snakemake@params[['index_threshold']]

dat <- fread(snakemake@input[[1]], sep = '\t', select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col))

dat[, (chr_col) := as.character(get(chr_col))]
dat[, (bp_col) := as.integer(get(bp_col))]

#dat <- dat[!(get(snp_col) %in% snakemake@params$snps_to_ignore)]

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
