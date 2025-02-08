library(data.table)

setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
p_col <- snakemake@config$p_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col

gwas <- fread(snakemake@input$gwas_file, sep = '\t', header = T, select = c(chr_col, bp_col, ref_col, alt_col, p_col, beta_col, se_col))
bim_ids <- fread(snakemake@input$range_file, sep = ' ', header = F)
names(bim_ids) <- 'SNPID'

gwas[, SNPID := paste(chr, bp, ref, alt, sep = ':'), env = list(chr = chr_col, bp = bp_col, ref = ref_col, alt = alt_col)]
gwas <- merge(gwas, bim_ids, by = 'SNPID')

# The SumHer documentation states that A1 is the 'test allele' and A2 is the 'other allele', so I take this to mean A1 is the effect allele and A2 the reference/other allele

setnames(gwas, c(chr_col, bp_col, alt_col, ref_col), c('chr', 'bp', 'A1', 'A2'))

#gwas <- gwas[!(chr %in% c('X', 'Y', 'MT'))]

gwas[, `:=` (len.A1 = nchar(A1), len.A2 = nchar(A2))]

gwas <- gwas[len.A1 == 1 & len.A2 == 1]

# 'Predictor' has format chr:bp:A2:A1 in my simgwas file, not sure it is correct, though
gwas[, Predictor := paste(chr, bp, A2, A1, sep = ':')]

gwas <- unique(gwas, by = 'Predictor')

gwas[, Z := beta/se, env = list(beta = beta_col, se = se_col)]
# qnorm(1e-100) = -21.27345
gwas[Z > 21, Z := 21]
gwas[Z < -21, Z := -21]
gwas <- na.omit(gwas, cols = c('Predictor', 'A1', 'A2', 'Z'))

gwas[, n := snakemake@params$N]

fwrite(gwas[, .(Predictor, A1, A2, n, Z)], sep = '\t', file =  snakemake@output[[1]])
