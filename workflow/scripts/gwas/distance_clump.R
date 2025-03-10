library(data.table)
setDTthreads(snakemake@threads)
library(magrittr)

#save.image('clump.RData')
#stop()

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
p_col <- snakemake@config$p_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
rsid_col <- snakemake@config$rsid_col

distance_window <- snakemake@params[['distance_window']]
index_threshold <- snakemake@params[['index_threshold']]

dat <- fread(snakemake@input[[1]], sep = '\t', select = c(chr_col, bp_col, ref_col, alt_col, rsid_col, p_col, beta_col, se_col))

if(!snakemake@params$mhc) {
  dat <- dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

dat[, `:=` (chr = as.character(chr), pos = as.integer(pos), p = as.numeric(p)), env = list(chr = chr_col, pos = bp_col, p = p_col)]

dat <- dat[p <= index_threshold, env = list(p = p_col)]

setorderv(dat, p_col)

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
