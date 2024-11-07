library(ggcorrplot)
library(dplyr)
library(magrittr)
library(data.table)

metadat <- fread(snakemake@input[['metadata']], select = c('abbrv', 'pretty_name'))

dat <- fread(snakemake@input[['rg']])

dat <- dat[!(trait.A %in% c('spondylo', 'igan-sans-ic')) & !(trait.B %in% c('spondylo', 'igan-sans-ic'))]

traits <- unique(dat$trait.B)

dat[, .(trait.A, trait.B, rg.sr)] %>%
  rbind(., .[, .(trait.A = trait.B, trait.B = trait.A, rg.sr)]) %>%
  tidyr::pivot_wider(names_from = trait.A, values_from = rg.sr, values_fill = 1) %>%
  data.table %>%
  .[trait.B %in% traits] %>%
  .[order(match(trait.B, traits))] %>%
  setcolorder(c('trait.B', traits)) %>%
  .[, ..traits] %>%
  as.matrix %>%
  set_rownames(., traits) -> rg_mat

metadat <- metadat[abbrv %in% traits][order(match(abbrv, traits))]
rownames(rg_mat) <- metadat$pretty_name
colnames(rg_mat) <- metadat$pretty_name

png(filename = snakemake@output[[1]], width = 10, height = 10, res = 400, units = 'in', type = 'cairo')
phmap <- pheatmap::pheatmap(rg_mat, cellwidth = 20, cellheight = 20, breaks = seq(-1, 1, length.out = 100), cutree_rows = 3, cutree_cols = 3)
dev.off()
