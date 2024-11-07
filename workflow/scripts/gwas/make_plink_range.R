library(data.table)
options(warn = 2)

setDTthreads(snakemake@threads)

gwas_file <- snakemake@input[['gwas_file']]

bim_file <- snakemake@input[['bim_file']]

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
ref_col <- snakemake@params[['ref_col']]
alt_col <- snakemake@params[['alt_col']]
tmpdir <- snakemake@resources[['tmpdir']]

gwas_dat <- fread(gwas_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col), tmpdir = tmpdir)

gwas_dat[ , (chr_col) := as.character(get(chr_col))]
gwas_dat[ , (bp_col) := as.integer(get(bp_col))]

bim_dat <- fread(bim_file, sep = '\t', header = F, col.names = c('CHR38', 'ID', 'Cm', 'BP38', 'A1', 'A2'))

bim_dat[, 'CHR38' := as.character(CHR38)]
bim_dat[, 'BP38' := as.integer(BP38)]

bim_join <- merge(bim_dat, gwas_dat, by.x = c(chr_col, bp_col), by.y = c('CHR38', 'BP38'), sort = F)

# Make sure alleles match
bim_join <- bim_join[(get(ref_col) == A1 & get(alt_col) == A2) | (get(ref_col) == A2 & get(alt_col) == A1)]

if(snakemake@params$mhc == F) {
  bim_join <- bim_join[!(CHR38 == 6 & BP38 %between% c(24e6, 45e6))]
}

fwrite(bim_join[, .(ID)], file = snakemake@output[[1]], row.names = F, sep = ' ', col.names = F, quote = F)
