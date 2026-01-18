library(ggplot2)
library(data.table)
library(ensembldb)
library(patchwork)
library(locuszoomr)

edb <- EnsDb(snakemake@input$edb)
ighkl_dat <- fread(snakemake@input$sumstats)

chrom <- as.character(snakemake@config$loci$igh$chrom)
locus_range <- c(snakemake@config$loci$igh$start, snakemake@config$loci$igh$stop)

draw_ig_locus_scatter_plot <- function(study, isotype) {
  if (study == "meta") {
    p_value_label <- sprintf("p_value.%s_%s", isotype, study)
  } else {
    p_value_label <- sprintf("p_value.%s.%s", study, isotype)
  }

  loc <- locus(na.omit(ighkl_dat, p_value_label), ens_db = edb, chrom = "chromosome", pos = "base_pair_location", p = p_value_label, labs = "rsid", seqname = chrom, xrange = locus_range)

  pl <- gg_scatter(loc, showLD = F, min.segment.length = 0)

  pl
}

pls <- list(
  draw_ig_locus_scatter_plot("meta", "igg"),
  draw_ig_locus_scatter_plot("meta", "iga"),
  draw_ig_locus_scatter_plot("meta", "igm")
)

study <- "meta"
p_value_label <- "p_value.igg_meta"

loc <- locus(na.omit(ighkl_dat, p_value_label), ens_db = edb, chrom = "chromosome", pos = "base_pair_location", p = p_value_label, labs = "rsid", seqname = chrom, xrange = locus_range)

pls[[length(pls)+1]] <- gg_genetracks(loc, filter_gene_biotype = c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene"))

p <- wrap_plots(pls, ncol = 1, heights = c(rep(1, length(pls) - 1), 1)) +
  plot_annotation(tag_levels = 'A')

ggsave(p, width = 6, height = 2.5 * length(pls), file = snakemake@output[[1]])
