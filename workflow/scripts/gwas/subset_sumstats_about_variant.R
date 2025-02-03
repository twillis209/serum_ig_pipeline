library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config[['chr_col']]
bp_col <- snakemake@config[['bp_col']]
ref_col <- snakemake@config[['ref_col']]
alt_col <- snakemake@config[['alt_col']]
p_col <- snakemake@config[['p_col']]
rsid_col <- snakemake@config[['rsid_col']]
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col
window <- as.numeric(snakemake@params[['window']])
variant_id <- snakemake@wildcards[['variant_id']]

cols <- c(chr_col, bp_col, ref_col, alt_col, beta_col, se_col, p_col, rsid_col)

dat <- fread(snakemake@input[[1]], sep = '\t', select = cols)

if(dat[get(rsid_col) == variant_id, .N] == 0) {
  stop("Couldn't find variant in summary statistics")
} else {
  variant_chr <- dat[get(rsid_col) == variant_id, get(chr_col)]
  variant_bp <- dat[get(rsid_col) == variant_id, get(bp_col)]

  dat <- dat[get(chr_col) == variant_chr]
  dat <- dat[get(bp_col) %between% c(variant_bp-window/2, variant_bp+window/2)]
  fwrite(dat[, ..cols], file = snakemake@output[['sum_stats']], sep = '\t')

  dat[, id_col := paste(chr_col, bp_col, ref_col, alt_col, sep = ':'), env = list(chr_col = chr_col, bp_col = bp_col, ref_col = ref_col, alt_col = alt_col)]

  fwrite(dat[, .(id_col)], file = snakemake@output[['ids']], sep = '\t', col.names = F)
}
