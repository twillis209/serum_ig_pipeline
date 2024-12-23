library(data.table)
options(warn = 2)

setDTthreads(snakemake@threads)

gwas_file <- snakemake@input[['gwas_file']]
metadata_file <- snakemake@input[['metadata_file']]
trait_A <- snakemake@params[['trait_A']]
trait_B <- snakemake@params[['trait_B']]
chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
beta_a_col <- snakemake@config$beta_a_col
beta_b_col <- snakemake@config$beta_b_col
se_a_col <- snakemake@config$se_a_col
se_b_col <- snakemake@config$se_b_col
gwas_file_A <- snakemake@output[['gwas_file_A']]
gwas_file_B <- snakemake@output[['gwas_file_B']]
tmpdir <- snakemake@resources[['tmpdir']]

gwas_dat <- fread(gwas_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, beta_a_col, beta_b_col, se_a_col, se_b_col), tmpdir = tmpdir)

meta_dat <- fread(metadata_file, sep = '\t', header = T)
trait_A_N <- meta_dat[abbrv == trait_A, N]
trait_A_N <- gsub(',', '', trait_A_N)

trait_B_N <- meta_dat[abbrv == trait_B, N]
trait_B_N <- gsub(',', '', trait_B_N)
# The SumHer documentation states that A1 is the 'test allele' and A2 is the 'other allele', so I take this to mean A1 is the minor allele and A2 the major allele

setnames(gwas_dat, c(chr_col, bp_col, alt_col, ref_col), c('chr', 'bp', 'A1', 'A2'))

gwas_dat <- gwas_dat[!(chr %in% c('X', 'Y', 'MT'))]

gwas_dat[, `:=` (len.A1 = nchar(A1), len.A2 = nchar(A2))]

gwas_dat <- gwas_dat[len.A1 == 1 & len.A2 == 1]

# 'Predictor' has format chr:bp:A2:A1 in my simgwas file, not sure it is correct, though
gwas_dat[, Predictor := paste(chr, bp, A2, A1, sep = ':')]

gwas_dat <- unique(gwas_dat, by = 'Predictor')

gwas_dat[, Z := get(beta_a_col)/get(se_a_col)]
# qnorm(1e-100) = -21.27345
gwas_dat[Z > 21, Z := 21]
gwas_dat[Z < -21, Z := -21]
gwas_dat[, n := trait_A_N]
gwas_dat <- na.omit(gwas_dat, cols = c('Predictor', 'A1', 'A2', 'Z'))

fwrite(gwas_dat[, .(Predictor, A1, A2, n, Z)], sep = '\t', file =  gwas_file_A)

gwas_dat[, Z := get(beta_b_col)/get(se_b_col)]
gwas_dat[Z > 21, Z := 21]
gwas_dat[Z < -21, Z := -21]
gwas_dat[, n := trait_B_N]
gwas_dat <- na.omit(gwas_dat, cols = c('Predictor', 'A1', 'A2', 'Z'))

fwrite(gwas_dat[, .(Predictor, A1, A2, n, Z)], sep = '\t', file =  gwas_file_B)
