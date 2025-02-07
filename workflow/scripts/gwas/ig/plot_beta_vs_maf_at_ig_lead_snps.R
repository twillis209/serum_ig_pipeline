library(data.table)
library(ggplot2)
library(magrittr)

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

cols <- c('rsid', 'beta', 'standard_error', 'gnomadNFE', 'Data set')

merged <- rbindlist(list(igm[, ..cols],
               igg[, ..cols],
               iga[, ..cols]),
               fill = T
               )

merged[, MAF := ifelse(gnomadNFE > 0.5, 1-gnomadNFE, gnomadNFE)]

pl <- ggplot(data = merged[!is.na(MAF)][MAF > 0])+
  geom_point(aes(x = MAF, y = abs(beta), col = `Data set`), alpha = 0.4)+
  xlab('MAF')+
  ylab('|Effect estimate|')+
  scale_x_log10(limits = c(1e-5, 1), breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))+
  coord_fixed()

ggsave(pl, file = snakemake@output[[1]], width = 4, height = 1.5)
