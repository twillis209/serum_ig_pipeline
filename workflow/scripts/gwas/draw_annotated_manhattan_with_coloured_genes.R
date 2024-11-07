library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(ggnewscale)
library(tidyverse)
library(pidPipelineCode)

theme_set(theme_bw())

colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
red <- "#009E73"
green <- "#D55E00"

manh_plot_width <- 370
manh_plot_height <- 210
gene_text_size <- 7
point_size <- 0.6

chrom_col <- snakemake@params[['chrom_col']]
bp_col <- snakemake@params[['bp_col']]
prin_col <- snakemake@params[['prin_col']]
snp_col <- snakemake@params[['snp_col']]
ylim <- ifelse(is.null(snakemake@params[['ylim']]), 20, snakemake@params[['ylim']])

gwas_dat <- fread(snakemake@input[['gwas']], sep = '\t', select = c(chrom_col, bp_col, prin_col, snp_col))

setnames(gwas_dat, c(chrom_col, bp_col, prin_col), c('chr', 'bp', 'p'))

gwas_dat <- gwas_dat[chr %in% seq(1, 22)]

gwas_dat[, chr := as.integer(chr)]

gwas_dat <- gwas_dat[!is.na(chr)]

munged_gwas_input <- munge_gwas_input(gwas_dat)

munged_gwas_input$gwas_data[, chr := as.character(chr)]

rsID_dat <- fread(snakemake@input[['rsIDs']], sep = '\t', select = c('SNPID', 'rsID', 'topGene', prin_col))
setnames(rsID_dat, prin_col, 'p')

rsID_dat[, c('chr', 'bp') := tstrsplit(SNPID, split = ':', keep = 1:2)]
rsID_dat[, `:=` (chr = as.character(chr), bp = as.integer(bp))]

merged_dat <- rbindlist(
  list(
munged_gwas_input$gwas_data[p < 5e-8][, .(chr, bp, bp_add, bp_cum, gwas.p = p)][rsID_dat[p < 5e-8], on = .(chr, bp), roll = 'nearest', nomatch = 0],
munged_gwas_input$gwas_data[p < 1e-5 & p > 5e-8][, .(chr, bp, bp_add, bp_cum, gwas.p = p)][rsID_dat[p < 1e-5 & p > 5e-8], on = .(chr, bp), roll = 'nearest', nomatch = 0]
  )
)

merged_dat[, p := NULL]
setnames(merged_dat, 'gwas.p', 'p')
merged_dat <- merged_dat[!(chr == '6' & bp %between% c(24e6, 45e6))]

#black_dat <- merged_dat[!(topGene %in% c(snakemake@params[['red_highlights']], snakemake@params[['green_highlights']]))]
red_dat <- merged_dat[topGene %in% snakemake@params[['red_highlights']]]
red_dat[, colour := red]
green_dat <- merged_dat[topGene %in% snakemake@params[['green_highlights']]]
green_dat[, colour := green]

if(red_dat[, .N] == 0 & green_dat[, .N] == 0) {
  repel_dat <- merged_dat
  repel_dat[, colour := 'black']
} else {
  repel_dat <- rbind(red_dat, green_dat)
}

manhplot <- ggplot(munged_gwas_input$gwas_data, aes(x = bp_cum, y = -log10(p), color = as_factor(chr), size = -log10(p))) +
  geom_hline(yintercept = -log10(5e-8), color = "grey40", linetype = "dashed") + 
  geom_hline(yintercept = -log10(1e-5), color = "grey40", linetype = "dashed") + 
  geom_point(size = 0.3) +
  scale_x_continuous(label = munged_gwas_input$axis_set$chr, breaks = munged_gwas_input$axis_set$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_size_continuous(range = c(0.5,3)) +
  scale_color_manual(values = rep(colorblind_palette, unique(length(munged_gwas_input$axis_set$chr))))+
  labs(x = NULL, y = "-log<sub>10</sub>(p)")+
  new_scale_color()+
  geom_text_repel(aes(label = topGene, col = colour), hjust = -0.2, data = repel_dat, size = gene_text_size, force = 1.5)+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.background = element_rect(color = 'transparent', fill = 'transparent'),
    axis.text.x=element_text(size= 16, color="black"),
    axis.text.y=element_text(size= 16, color="black"),
    legend.position = 'none',
    axis.title.y = element_markdown(angle = 0, size = 24, vjust = 0.5))

ggsave(manhplot, file = snakemake@output[[1]], width = manh_plot_width, height = manh_plot_height, unit = 'mm')
