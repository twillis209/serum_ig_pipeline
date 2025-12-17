library(ggplot2)
library(data.table)
library(ensembldb)
library(patchwork)
library(locuszoomr)

edb <- EnsDb(snakemake@input$edb)
ighkl_dat <- fread(snakemake@input$sumstats)

draw_ig_locus_scatter_plot <- function(study, isotype, ig_locus, with_gene_tracks = FALSE) {
  p_value_label <- sprintf("p_value.%s.%s", study, isotype)

  loc <- locus(na.omit(ighkl_dat, p_value_label), ens_db = edb, chrom = "chromosome", pos = "base_pair_location", p = p_value_label, labs = "rsid", seqname = as.character(coords[[ig_locus]]$chrom), xrange = c(coords[[ig_locus]]$start, coords[[ig_locus]]$stop))

  pl <- gg_scatter(loc, showLD = F, min.segment.length = 0) + ggtitle(sprintf("%s.%s", study, isotype))

  if (with_gene_tracks) pl <- pl / gg_genetracks(loc)

  pl
}

lapply(snakemake@config[[paste0(snakemake@wildcards$isotype, "_studies")]], draw_ig_locus_scatter_plot, isotype = snakemake@wildcards$isotype, ig_locus = snakemake@wildcards$ighkl_locus)



# Reduce? applying the patchwork `/`?

ggsave(file = snakemake@output$igh)
ggsave(file = snakemake@output$igk)
ggsave(file = snakemake@output$igl)
