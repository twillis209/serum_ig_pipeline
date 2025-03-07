library(data.table)
library(ggplot2)
library(patchwork)
library(magrittr)
library(serumIgPipelineCode)

theme_set(theme_bw())

dat <- fread(snakemake@input[[1]])

dat[, `:=` (z1 = b1/se1, z2 = b2/se2),
    env = list(b1 = sprintf("beta.%s", snakemake@wildcards$first_isotype),
               se1 = sprintf("standard_error.%s", snakemake@wildcards$first_isotype),
               b2 = sprintf("beta.%s", snakemake@wildcards$second_isotype),
               se2 = sprintf("standard_error.%s", snakemake@wildcards$second_isotype)
               )
   ]

melt(dat, measure.vars = measure(stat, dataset, sep = '.'))[stat == 'p_value'] %>%
  ggplot()+
  geom_point(aes(x = base_pair_location, y = value))+
  scale_y_neglog10()+
  ylab('-log10(p)')+
  facet_grid(dataset ~ .) -> locus_plot

zscore_plot <- ggplot(dat)+
  geom_point(aes(x = z1, y = z2))+
  xlab(snakemake@wildcards$first_isotype)+
  ylab(snakemake@wildcards$second_isotype)+
  geom_vline(xintercept = qnorm(2.5e-8), linetype = 'dashed', col = 'blue')+
  geom_vline(xintercept = qnorm(2.5e-8, lower.tail = F), linetype = 'dashed', col = 'blue')+
  geom_hline(yintercept = qnorm(2.5e-8), linetype = 'dashed', col = 'blue')+
  geom_hline(yintercept = qnorm(2.5e-8, lower.tail = F), linetype = 'dashed', col = 'blue')+
  coord_fixed(ratio = 1)

ggsave(locus_plot+zscore_plot+plot_annotation(title = sprintf('%s and %s', snakemake@wildcards$first_rsid, snakemake@wildcards$second_rsid)), width = 6, height = 4, file = snakemake@output[[1]])
