library(data.table)
setDTthreads(snakemake@threads)

save.image('merge_r2_with_sumstats.RData')
stop()

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
rsid_col <- snakemake@config$rsid_col

gwas <- fread(snakemake@input[['gwas']])

coord_cols <- c(chr_col, bp_col, ref_col, alt_col)

variant_coord_dat <- gwas[get(rsid_col) == snakemake@wildcards$variant_id, ..coord_cols]
variant_coord_dat[, id := paste(get(chr_col), get(bp_col), get(ref_col), get(alt_col), sep = ':')]
variant_coord_dat[, chr := as.character(chr), env = list(chr = chr_col)]
variant_coord_dat[, bp := as.integer(bp), env = list(bp = bp_col)]

ld <- fread(snakemake@input[['ld']], header = F)

ld_vars <- fread(snakemake@input[['ld_vars']], header = F, col.names = 'ID')

ld_vars[, c('chr', 'pos', 'ref', 'alt') := tstrsplit(ID, split = ':')]
ld_vars[, chr := as.character(chr)]
ld_vars[, pos := as.integer(pos)]

merged <- merge(variant_coord_dat, ld_vars, by.x = c(chr_col, bp_col), by.y = c('chr', 'pos'))

if(merged[, .N] == 0) {
  stop("Lead SNP not present in LD matrix")
} else if(merged[!((ref == get(ref_col) & alt == get(alt_col)) | (ref == get(alt_col) & alt == get(ref_col)))]) {
  stop("Lead SNP not present in LD matrix")
} else {
  names(ld) <- ld_vars$ID

  ld[, ID := ld_vars$ID]

  molten_ld <- melt(ld[ID == variant_coord_dat[, id]], id.vars = 'ID')
  molten_ld[, ID := NULL]
  names(molten_ld) <- c('SNPID', 'r2')

  merged <- merge(gwas, molten_ld, by = 'SNPID', all.x = T)

  fwrite(merged, file = snakemake@output[[1]], sep = '\t')
}
