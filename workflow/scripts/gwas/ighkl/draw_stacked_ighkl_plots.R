library(ggplot2)
library(data.table)
library(ensembldb)
library(patchwork)
library(locuszoomr)

edb <- EnsDb(snakemake@input$edb)
ighkl_dat <- fread(snakemake@input$sumstats)

draw_ig_locus_scatter_plot <- function(study, isotype, ig_locus) {
  p_value_label <- sprintf("p_value.%s.%s", study, isotype)

  loc <- locus(na.omit(ighkl_dat, p_value_label), ens_db = edb, chrom = "chromosome", pos = "base_pair_location", p = p_value_label, labs = "rsid", seqname = snakemake@config$loci[[ig_locus]]$chrom, xrange = c(snakemake@config$loci[[ig_locus]]$start, snakemake@config$loci[[ig_locus]]$stop))

  pl <- gg_scatter(loc, showLD = F, min.segment.length = 0) + ggtitle(sprintf("%s, %s", snakemake@config$gwas_datasets[[paste(study, isotype, sep = "-")]]$pretty_study, snakemake@config$pretty_isotypes[[isotype]]))

  pl
}

pls <- lapply(snakemake@config[[paste0(snakemake@wildcards$isotype, "_studies")]], draw_ig_locus_scatter_plot, isotype = snakemake@wildcards$isotype, ig_locus = snakemake@wildcards$ighkl_locus)

study <- snakemake@config[[paste0(snakemake@wildcards$isotype, "_studies")]][1]
p_value_label <- sprintf("p_value.%s.%s", study, snakemake@wildcards$isotype)
ig_locus <- snakemake@wildcards$ighkl_locus

loc <- locus(na.omit(ighkl_dat, p_value_label), ens_db = edb, chrom = "chromosome", pos = "base_pair_location", p = p_value_label, labs = "rsid", seqname = snakemake@config$loci[[ig_locus]]$chrom, xrange = c(snakemake@config$loci[[ig_locus]]$start, snakemake@config$loci[[ig_locus]]$stop))

pls[[length(pls)+1]] <- gg_genetracks(loc, filter_gene_biotype = c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene"))

p <- wrap_plots(pls, ncol = 1, heights = c(rep(1, length(pls)-1), 3))

ggsave(p, width = 6, height = 2.5 * length(pls), file = snakemake@output[[1]])
