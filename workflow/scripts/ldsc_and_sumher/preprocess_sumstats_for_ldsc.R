library(data.table)
setDTthreads(snakemake@threads)

save.image('preprocess_sumstats.RData')

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
p_col <- snakemake@config$p_col
beta_col <- snakemake@config$beta_col
snp_col <- snakemake@config$id_col

sumstats <- fread(snakemake@input[['sumstats']], sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, snp_col))
sumstats[, (chr_col) := as.character(get(chr_col))]
maf <- fread(snakemake@input[['maf']], select = c('ID', 'ALT_FREQS'))

maf[, c('chr', 'bp', 'ref', 'alt') := tstrsplit(ID, split = ':')]
maf[, chr := as.character(chr)]
maf[, bp := as.integer(bp)]

merged <- merge(sumstats, maf, by.y = c('chr', 'bp', 'ref', 'alt'), by.x = c(chr_col, bp_col, ref_col, alt_col))

# Need SNP column which can be merged with that in the LD scores file
merged[, (snp_col) := paste(get(chr_col), get(bp_col), get(ref_col), get(alt_col), sep = ':')]

cols <- c(snp_col, chr_col, bp_col, ref_col, alt_col, p_col, beta_col, 'ALT_FREQS')

fwrite(merged[, ..cols], sep = '\t', file = snakemake@output[[1]])
