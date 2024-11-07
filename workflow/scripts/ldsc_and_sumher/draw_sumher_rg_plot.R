library(data.table)
library(ggplot2)
library(patchwork)
library(pidPipelineCode)

theme_set(theme_bw()+
          theme(
            axis.text.x = element_text(angle = 90)
          ))

rg <- fread(snakemake@input[['rg']])

metadat <- fread(snakemake@input[['metadata']])

gps <- fread(snakemake@input[['gps']])

gps <- gps[!(trait_B %in% c('spondylo', 'igan-sans-ic'))]

rg <- rg[!(trait.B %in% c('spondylo', 'igan-sans-ic'))]

merged <- merge(rg, metadat[, .(abbrv, pretty_name)], by.x = 'trait.B', by.y = 'abbrv', all.x = T)

merged <- merge(merged, gps[, .(trait_B, gps.p = pvalue)], by.x = 'trait.B', by.y = 'trait_B')

merged[, `:=` (rg.fdr = p.adjust(rg.p.sr, method = 'BH'), gps.fdr = p.adjust(gps.p, method = 'BH'))]

rg_pl <- ggplot(merged)+
  geom_point(aes(y = reorder(pretty_name, rg.p.sr), x = rg.sr, col = rg.fdr <= 0.05))+
  geom_errorbarh(aes(y = reorder(pretty_name, rg.p.sr), xmin = rg.sr-1.96*rg.se.sr, xmax = rg.sr+1.96*rg.se.sr, col = rg.fdr <= 0.05))+
  geom_vline(xintercept = 0, linetype = 'dashed')+
  ylab('Trait')+
  xlab('rg')+
  theme(legend.position = 'none')

gps_pl <- ggplot(merged)+
  geom_point(aes(y = reorder(pretty_name, rg.p.sr), x = gps.p, col = gps.fdr <= 0.05))+
  ylab('Trait')+
  xlab('GPS p-value')+
  scale_x_neglog10()+
  theme(legend.position = 'none')


ggsave(rg_pl+gps_pl, file = snakemake@output[[1]], width = 12, height = 10)
