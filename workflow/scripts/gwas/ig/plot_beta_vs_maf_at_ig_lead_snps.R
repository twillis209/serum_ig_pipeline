library(data.table)
library(ggplot2)
library(magrittr)
library(patchwork)
library(scales)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=6),
            plot.title=element_text(hjust=0.5, size=10),
            strip.text=element_text(size=6),
            axis.text.x=element_text(size=6, color="black"),
            axis.text.y=element_text(size=6, color="black"),
            legend.title=element_text(size=6),
            legend.text=element_text(size=6)
          )
)

igm <- fread(snakemake@input$igm)
igm[, `Data set` := 'IgM']
igg <- fread(snakemake@input$igg)
igg[, `Data set` := 'IgG']
iga <- fread(snakemake@input$iga)
iga[, `Data set` := 'IgA']

cols <- c('rsid', 'beta', 'standard_error', 'maf.meta', 'Data set')

merged <- rbindlist(list(igm[, ..cols],
               igg[, ..cols],
               iga[, ..cols]),
               fill = T
               )

pl1 <- ggplot(data = merged[!is.na(maf.meta)][maf.meta > 0])+
  geom_point(aes(x = maf.meta, y = abs(beta), col = `Data set`), alpha = 0.4)+
  xlab('MAF')+
  ylab('|Effect estimate|')+
  scale_x_continuous(limits = c(1e-2, 0.5), breaks = c(1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5))+
  ylim(c(0,1.2))


pl2 <- ggplot(data = merged[!is.na(maf.meta)][maf.meta > 0])+
  geom_point(aes(x = maf.meta, y = abs(beta), col = `Data set`), alpha = 0.4)+
  xlab('MAF')+
  ylab('|Effect estimate|')+
  scale_x_log10(limits = c(1e-4, 1e-2), breaks = c(1e-4, 1e-3, 1e-2), labels = function(x) format(x, scientific = F))+
  coord_fixed()+
  theme(legend.position = 'none')

ggsave(pl2+pl1+plot_layout(axes = 'collect')+plot_annotation(tag_levels = 'A'), file = snakemake@output[[1]], width = 4, height = 1.5)
