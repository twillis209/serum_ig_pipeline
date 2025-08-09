library(data.table)
library(locuszoomr)
library(ensembldb)
library(ggplot2)
library(patchwork)

edb <- EnsDb(snakemake@input$edb)

dat <- fread(snakemake@input[[1]], sep = '\t', header = T, select = c('rsid', 'chromosome', 'base_pair_location', 'p_value'))

dat <- dat[chromosome == snakemake@wildcards$chrom & base_pair_location %between% c(snakemake@wildcards$start, snakemake@wildcards$end)]

dat <- na.omit(dat, c('p_value', 'base_pair_location'))

loc <- locus(data = dat, chrom = 'chromosome', pos = 'base_pair_location', p = 'p_value', labs = 'rsid', seqname = snakemake@wildcards$chrom, xrange = c(as.numeric(snakemake@wildcards$start), as.numeric(snakemake@wildcards$end)), ens_db = edb)

ggsave(gg_scatter(loc, showLD = F)/gg_genetracks(loc), file = snakemake@output[[1]], width = 6, height = 7)
