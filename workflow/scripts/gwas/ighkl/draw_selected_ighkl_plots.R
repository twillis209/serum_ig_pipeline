library(ggplot2)
library(data.table)
library(ensembldb)
library(patchwork)
library(locuszoomr)

edb <- EnsDb(snakemake@input$edb)
ighkl_dat <- fread(snakemake@input$sumstats)

get_loc <- function(study, isotype, ig_locus) {
  p_value_label <- sprintf("p_value.%s.%s", study, isotype)

  loc <- locus(na.omit(ighkl_dat, p_value_label), ens_db = edb, chrom = "chromosome", pos = "base_pair_location", p = p_value_label, labs = "rsid", seqname = snakemake@config$loci[[ig_locus]]$chrom, xrange = c(snakemake@config$loci[[ig_locus]]$start, snakemake@config$loci[[ig_locus]]$stop))
  loc
}

draw_ig_locus_scatter_plot <- function(study, isotype, ig_locus) {
  p_value_label <- sprintf("p_value.%s.%s", study, isotype)

  loc <- get_loc(study, isotype, ig_locus)

  pl <- gg_scatter(loc, showLD = F, min.segment.length = 0) + ggtitle(sprintf("%s, %s", snakemake@config$gwas_datasets[[paste(study, isotype, sep = "-")]]$pretty_study, snakemake@config$pretty_isotypes[[isotype]]))

  pl
}

epic_iga_igh_loc <- get_loc("epic", "iga", "igh")

iga_igh_pl <- wrap_plots(list(draw_ig_locus_scatter_plot("epic", "iga", "igh"),
  draw_ig_locus_scatter_plot("pietzner", "iga", "igh"),
  draw_ig_locus_scatter_plot("eldjarn", "iga", "igh"),
  gg_genetracks(epic_iga_igh_loc, filter_gene_biotype = c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene"))),
  heights = c(1, 1, 1, 3),
  ncol = 1
)

igg_igh_pl <- wrap_plots(list(draw_ig_locus_scatter_plot("epic", "igg", "igh"),
  draw_ig_locus_scatter_plot("pietzner", "igg", "igh"),
  draw_ig_locus_scatter_plot("eldjarn", "igg", "igh"),
  gg_genetracks(epic_iga_igh_loc, filter_gene_biotype = c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene"))),
  heights = c(1, 1, 1, 3),
  ncol = 1
)

igm_igh_pl <- wrap_plots(list(
  draw_ig_locus_scatter_plot("eldjarn", "igg", "igh"),
  gg_genetracks(epic_iga_igh_loc, filter_gene_biotype = c("protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene"))),
  heights = c(1, 3),
  ncol = 1
)
