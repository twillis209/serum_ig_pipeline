library(data.table)

gwas_file_a <- snakemake@input[['A']]
gwas_file_b <- snakemake@input[['B']]
chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
p_col <- snakemake@config$p_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
id_col <- snakemake@config$id_col
output_file <- snakemake@output[['AB']]
mhc <- snakemake@params[['mhc']]
join <- snakemake@params[['join']]
tmpdir <- snakemake@resources[['tmpdir']]

setDTthreads(snakemake@threads)

dat_a <- fread(gwas_file_a, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col, id_col), tmpdir = tmpdir)
dat_b <- fread(gwas_file_b, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col, id_col), tmpdir = tmpdir)

dat_a[ , c(ref_col, alt_col) := list(toupper(get(ref_col)), toupper(get(alt_col)))]
dat_a <- dat_a[get(ref_col) %in% c('A','T','C','G') & get(alt_col) %in% c('A','T','C','G')]
dat_a[, (chr_col) := as.character(get(chr_col))]
dat_a <- na.omit(dat_a)

dat_b[ , c(ref_col, alt_col) := list(toupper(get(ref_col)), toupper(get(alt_col)))]
dat_b <- dat_b[get(ref_col) %in% c('A','T','C','G') & get(alt_col) %in% c('A','T','C','G')]
dat_b[, (chr_col) := as.character(get(chr_col))]
dat_b <- na.omit(dat_b)

if(join == 'inner') {
  merged_dat <- merge(dat_a, dat_b, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else if(join == 'left') {
  merged_dat <- merge(dat_a, dat_b, all.x = T, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else if(join == 'right') {
  merged_dat <- merge(dat_a, dat_b, all.y = T, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else if(join == 'outer') {
  merged_dat <- merge(dat_a, dat_b, all.x = T, all.y = T, by = c(chr_col, bp_col), suffixes = c('.A', '.B'))
} else {
  stop(sprintf("Unrecognised join param: %s", join))
}

# Removes the MHC
if(!mhc) {
  merged_dat <- merged_dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
}

ref_a <- paste0(ref_col, '.A')
ref_b <- paste0(ref_col, '.B')
alt_a <- paste0(alt_col, '.A')
alt_b <- paste0(alt_col, '.B')

# Handle flipped alleles
                                        # TODO need to handle case where both not present
if(join == 'inner') {
  merged_dat <- merged_dat[
  (get(ref_a) == get(ref_b) & get(alt_a) == get(alt_b)) |
  (get(ref_a) == get(alt_b) & get(alt_a) == get(ref_b))
]
} else if(join == 'outer') {
  merged_dat <- merged_dat[
  (get(ref_a) == get(ref_b) & get(alt_a) == get(alt_b)) |
  (get(ref_a) == get(alt_b) & get(alt_a) == get(ref_b)) |
  (is.na(get(ref_a)) & is.na(get(alt_a))) |
  (is.na(get(ref_b)) & is.na(get(alt_b)))
  ]
} else {
  stop("Not yet implemented for this join type")
}

beta_to_change <- paste0(beta_col, '.B')

merged_dat[(get(ref_a) == get(alt_b) & get(alt_a) == get(ref_b)), (beta_to_change) := -get(beta_to_change)]

merged_dat[is.na(REF.A) & is.na(ALT.A), `:=` (REF.A = REF.B, ALT.A = ALT.B)]
merged_dat[is.na(SNPID.A), SNPID.A := SNPID.B]

merged_dat[, c(ref_b, alt_b) := NULL]
setnames(merged_dat, c(ref_a, alt_a), c(ref_col, alt_col))

merged_dat <- unique(merged_dat, by = c(chr_col, bp_col))

fwrite(merged_dat, file = output_file, sep = '\t', row.names = F)
