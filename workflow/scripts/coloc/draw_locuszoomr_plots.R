library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(gridExtra)

theme_set(theme_bw())

dat <- fread(snakemake@input[[1]])

seqname <- dat[, unique(chromosome)]
min_pos <- dat[, min(base_pair_location)]
max_pos <- dat[, max(base_pair_location)]

if(!is.null(snakemake@wildcards$first_isotype)) {
  first_trait_label <- snakemake@wildcards$first_isotype
  second_trait_label <- snakemake@wildcards$second_isotype
} else if(!is.null(snakemake@wildcards$non_ig)) {
  first_trait_label <- snakemake@wildcards$isotype
  second_trait_label <- snakemake@wildcards$non_ig
} else {
  stop("Can't find wildcard in input file path to identify summary statistics' suffix labels.")
}

dat <- na.omit(dat, c(sprintf('p_value.%s', first_trait_label), sprintf('p_value.%s', second_trait_label)))

dat[, `:=` (z1 = b1/se1, z2 = b2/se2),
    env = list(b1 = sprintf("beta.%s", first_trait_label),
               se1 = sprintf("standard_error.%s", first_trait_label),
               b2 = sprintf("beta.%s", second_trait_label),
               se2 = sprintf("standard_error.%s", second_trait_label)
               )
    ]

z_min <- dat[, min(z1, z2)]
z_max <- dat[, max(z1, z2)]

loc_A <- locus(data = dat, chrom = 'chromosome', pos = 'base_pair_location', p = sprintf('p_value.%s', first_trait_label), labs = 'rsid', seqname = seqname, xrange = c(min_pos, max_pos), ens_db = 'EnsDb.Hsapiens.v86')
loc_B <- locus(data = dat, chrom = 'chromosome', pos = 'base_pair_location', p = sprintf('p_value.%s', second_trait_label), labs = 'rsid', seqname = seqname, xrange = c(min_pos, max_pos), ens_db = 'EnsDb.Hsapiens.v86')

pls <- list()

pls[[1]] <- gg_scatter(loc_A, labels = snakemake@wildcards$first_rsid, showLD = F, min.segment.length = 0, nudge_y = 5)
pls[[2]] <- gg_scatter(loc_B, labels = snakemake@wildcards$second_rsid, showLD = F, min.segment.length = 0, nudge_y = 5)
pls[[3]] <- gg_genetracks(loc_A)
pls[[4]] <- ggplot(dat)+
  geom_point(aes(x = z1, y = z2))+
  xlab(first_trait_label)+
  ylab(second_trait_label)+
  geom_vline(xintercept = qnorm(2.5e-8), linetype = 'dashed', col = 'blue')+
  geom_vline(xintercept = qnorm(2.5e-8, lower.tail = F), linetype = 'dashed', col = 'blue')+
  geom_hline(yintercept = qnorm(2.5e-8), linetype = 'dashed', col = 'blue')+
  geom_hline(yintercept = qnorm(2.5e-8, lower.tail = F), linetype = 'dashed', col = 'blue')+
  coord_fixed(ratio = 1)+
  xlim(c(z_min, z_max))+
  ylim(c(z_min, z_max))

ggsave(arrangeGrob(grobs = pls, layout_matrix = cbind(1:3, 4),
                   top = sprintf('%s:%s-%s',
                                 snakemake@params$chrom,
                                 format(min_pos, big.mark = ',', scientific = F),
                                 format(max_pos, big.mark = ',', scientific = F)
                                 )
                   ), file = snakemake@output[[1]], width = 12, height = 9, unit = 'in')
