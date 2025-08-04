library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(tomics)
library(patchwork)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            axis.text.x=element_text(size=6, angle=90, color="black"),
            axis.text.y=element_text(size=10, angle = 90, color="black"),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10),
            panel.grid.major.y = element_blank()
          )
          )

iga <- fread(snakemake@input$iga, select = c('chromosome', 'base_pair_location', 'p_value'))
igg <- fread(snakemake@input$igg, select = c('chromosome', 'base_pair_location', 'p_value'))
igm <- fread(snakemake@input$igm, select = c('chromosome', 'base_pair_location', 'p_value'))

iga_procd_sumstats <- process_sumstats_for_manhattan(iga)
iga_procd_sumstats$axis_set$chr <- c(as.character(1:22), 'X')
rm(iga)
igg_procd_sumstats <- process_sumstats_for_manhattan(igg)
igg_procd_sumstats$axis_set$chr <- c(as.character(1:22), 'X')
rm(igg)
igm_procd_sumstats <- process_sumstats_for_manhattan(igm)
igm_procd_sumstats$axis_set$chr <- c(as.character(1:22), 'X')
rm(igm)

iga_plot <- draw_manhattan(iga_procd_sumstats, stat_col = "p_value", y_limits = snakemake@params$iga_ylim) +
  geom_hline(yintercept = 5e-8, linetype = 'dashed', col = 'black')

igg_plot <- draw_manhattan(igg_procd_sumstats, stat_col = "p_value", y_limits = snakemake@params$igg_ylim) +
  geom_hline(yintercept = 5e-8, linetype = 'dashed', col = 'black')


igm_plot <- draw_manhattan(igm_procd_sumstats, stat_col = "p_value", y_limits = snakemake@params$igm_ylim) +
  geom_hline(yintercept = 5e-8, linetype = 'dashed', col = 'black')

ggsave((iga_plot / igg_plot / igm_plot)+plot_annotation(tag_levels = 'A'), file = snakemake@output[[1]], width = snakemake@params$width, height = snakemake@params$height)
