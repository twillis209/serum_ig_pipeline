library(data.table)

# If tmp fills up, get truncated files which don't produce errors; that's bad!
options(warn = 2)

gwas_file_a <- snakemake@input[['A']]
gwas_file_b <- snakemake@input[['B']]
chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
p_col <- snakemake@config$p_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
rsid_col <- snakemake@config$rsid_col
join <- snakemake@params[['join']]

setDTthreads(snakemake@threads)

dat_a <- fread(gwas_file_a, sep = '\t', header = T, select = c(rsid_col, chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col))
dat_b <- fread(gwas_file_b, sep = '\t', header = T, select = c(rsid_col, chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col))

suffix_A <- sprintf('.%s', snakemake@wildcards$trait_A)
suffix_B <- sprintf('.%s', snakemake@wildcards$trait_B)

coord_cols <- c(chr_col, bp_col, ref_col, alt_col)

dat_a[, chr := as.character(chr), env = list(chr = chr_col)]
dat_b[, chr := as.character(chr), env = list(chr = chr_col)]

if(join == 'inner') {
  merged_dat <- merge(dat_a, dat_b, by = coord_cols, suffixes = c(suffix_A, suffix_B))
} else if(join == 'left') {
  merged_dat <- merge(dat_a, dat_b, all.x = T, by = coord_cols, suffixes = c(suffix_A, suffix_B))
} else if(join == 'right') {
  merged_dat <- merge(dat_a, dat_b, all.y = T, by = coord_cols, suffixes = c(suffix_A, suffix_B))
} else if(join == 'outer') {
  merged_dat <- merge(dat_a, dat_b, all.x = T, all.y = T, by = coord_cols, suffixes = c(suffix_A, suffix_B))
} else {
  stop(sprintf("Unrecognised join param: %s", join))
}

# Removes the MHC
if(!snakemake@params$mhc) {
  merged_dat <- merged_dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

if (!snakemake@params$ighkl) {
  for (x in snakemake@config$loci) {
    merged_dat <- merged_dat[!(
      chr == as.character(x$chrom) &
      pos %between% c(as.integer(x$start), as.integer(x$stop))
    ), env = list(chr = chr_col, pos = bp_col)]
  }
}

rsid_col_a <- paste0(rsid_col, suffix_A)
rsid_col_b <- paste0(rsid_col, suffix_B)

merged_dat[, rsid := fcoalesce(rsid_a, rsid_b), env = list(rsid_a = rsid_col_a, rsid_b = rsid_col_b)]
merged_dat[, c(rsid_col_a, rsid_col_b) := NULL]


fwrite(merged_dat, file = snakemake@output$AB, sep = '\t', row.names = F)
