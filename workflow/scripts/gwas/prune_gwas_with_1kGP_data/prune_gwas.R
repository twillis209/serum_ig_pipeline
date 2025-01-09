library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col

gwas_dat <- fread(snakemake@input[['gwas_file']], sep = '\t', header = T)

gwas_dat[, (chr_col) := as.character(get(chr_col))]

prune_dat <- fread(snakemake@input[['prune_file']], sep = ' ', header = F)

prune_dat[, c('chr', 'bp', 'REF', 'ALT') := tstrsplit(V1, ':')]

prune_dat[, chr := as.character(chr)]
prune_dat[, bp := as.integer(bp)]

maf_dat <- fread(snakemake@input$maf_file, select = c('ID', 'ALT_FREQS'))

prune_dat <- merge(prune_dat, maf_dat, by.x = 'V1', by.y = 'ID', all.x = T)

merged_dat <- merge(gwas_dat, prune_dat, by.x = c(chr_col, bp_col, ref_col, alt_col), by.y = c('chr', 'bp', 'REF', 'ALT'))

fwrite(merged_dat, file = snakemake@output[[1]], sep = '\t', col.names = T, row.names = F, quote = F)
