library(data.table)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(ggplot2)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
p_col <- snakemake@config$p_col
rsid_col <- snakemake@config$rsid_col

dat <- fread(snakemake@input[[1]], sep = '\t', header = T, select = c(rsid_col, chr_col, bp_col, p_col))

variant_chr <- dat[get(rsid_col) == snakemake@wildcards$variant_id, get(chr_col)]
variant_pos <- dat[get(rsid_col) == snakemake@wildcards$variant_id, get(bp_col)]

loc <- locus(data = dat, chrom = chr_col, pos = bp_col, p = p_col, labs = rsid_col, seqname = variant_chr, xrange = c(variant_pos - snakemake@params$window/2, variant_pos + snakemake@params$window/2), ens_db = 'EnsDb.Hsapiens.v86')

ggsave(gg_scatter(loc, showLD = F), file = snakemake@output[[1]], width = 6, height = 4)
