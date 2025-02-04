library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(serumIgPipelineCode)

if(is.null(snakemake@params[['width']])) {
  width <- 8
} else {

  width <- snakemake@params[['width']]
}

if(is.null(snakemake@params[['height']])) {
  height <- 2
} else {
  height <- snakemake@params[['height']]
}

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            axis.text.x=element_text(size=6, angle=90, color="black"),
            axis.text.y=element_text(size=10, color="black"),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10),
            panel.grid.major.y = element_blank()
          )
          )

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
p_col <- snakemake@config$p_col

gwas_dat <- fread(snakemake@input[['gwas']], sep = '\t', select = c(chr_col, bp_col, p_col))

procd_sumstats <- process_sumstats_for_manhattan(gwas_dat)


if(!is.null(snakemake@params$y_axis_break)) {
  manh_plot <- draw_manhattan(procd_sumstats, y_axis_break = as.numeric(snakemake@params$y_axis_break))
} else {
  manh_plot <- draw_manhattan(procd_sumstats)
}

if(!is.null(snakemake@params$ylim)) {
  manh_plot <- manh_plot + ylim(as.numeric(snakemake@params$ylim))
}

if(!is.null(snakemake@input$lead_snps)) {
  lead_snps <- fread(snakemake@input$lead_snps, header = T)

  if(lead_snps[, .N] > 0) {
    merged <- merge(lead_snps, procd_sumstats$data,
                        by.x = c(chr_col, bp_col),
                        by.y = c('chr', 'bp')
                        )

    manh_plot <- manh_plot+
      geom_point(size = 0.9, pch = 21, colour = 'black', data = merged)+
      geom_text_repel(aes(label = topGene), hjust = -0.2, size = 4, data = merged, colour = 'black')
    }
}

ggsave(manh_plot, file = snakemake@output[[1]], width = width, height = height)
