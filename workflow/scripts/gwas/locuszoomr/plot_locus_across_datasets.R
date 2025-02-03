library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(gridExtra)
library(patchwork)
library(stringr)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col

dat <- fread(snakemake@input[[1]], sep = '\t', header = T, select = c(chr_col, bp_col, snakemake@params$p_value_cols))

dat <- dat[chr == snakemake@params$chrom & pos %between% c(snakemake@params$start_pos, snakemake@params$stop_pos), env = list(chr = chr_col, pos = bp_col)]

dat[, snpid := paste(chr, bp), env = list(chr = chr_col, bp = bp_col)]

pls <- list()

pdf(snakemake@output[[1]], width = 7, height = 10)

for(i in seq_along(snakemake@params$p_value_cols)) {
  p_value_col <- snakemake@params$p_value_cols[[i]]

  loc <- locus(data = dat[!is.na(get(p_value_col))], chrom = chr_col, pos = bp_col, p = p_value_col, ens_db = "EnsDb.Hsapiens.v86", seqname = as.integer(snakemake@params$chrom), xrange = c(snakemake@params$start_pos, snakemake@params$stop_pos), labs = 'snpid')

  print(str_split_i(p_value_col, '\\.', 2))

  pl <- gg_scatter(loc, showLD = F, min.segment.length = 0, nudge_y = 5)/
    gg_genetracks(loc)+
    plot_annotation(title = str_split_i(p_value_col, '\\.', 2))

  print(pl)
}

dev.off()
