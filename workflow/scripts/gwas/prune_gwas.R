library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]

gwas_dat <- fread(snakemake@input[['gwas_file']], sep = '\t', header = T)

gwas_dat[, (chr_col) := as.character(get(chr_col))]

prune_dat <- fread(snakemake@input[['prune_file']], sep = ' ', header = F)

prune_dat[, c('chr', 'bp', 'REF', 'ALT') := tstrsplit(V1, ':')]

prune_dat[, chr := as.character(chr)]
prune_dat[, bp := as.integer(bp)]

merged_dat <- merge(gwas_dat, prune_dat, by.x = c(chr_col, bp_col, ref_col, alt_col), by.y = c('chr', 'bp', 'REF', 'ALT'))

fwrite(merged_dat, file = snakemake@output[[1]], sep = '\t', col.names = T, row.names = F, quote = F)
