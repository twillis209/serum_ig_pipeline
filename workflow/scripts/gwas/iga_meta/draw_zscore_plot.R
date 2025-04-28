library(data.table)
setDTthreads(snakemake@threads)
library(ggplot2)
library(stringr)

theme_set(theme_bw() +
        theme(
          axis.title = element_text(size=12),
          axis.text = element_text(size = 10, color = "black"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10)
        )
        )

ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
beta_a_col <- snakemake@params$beta_a
beta_b_col <- snakemake@params$beta_b
se_a_col <- snakemake@params$se_a
se_b_col <- snakemake@params$se_b

dat <- fread(snakemake@input[[1]], select = c('chromosome', 'base_pair_location', beta_a_col, beta_b_col, se_a_col, se_b_col))

dat <- na.omit(dat, c(beta_a_col, beta_b_col, se_a_col, se_b_col))

if(snakemake@wildcards$ighkl_inclusion == 'sans_ighkl') {
  for(x in c('igh', 'igk', 'igl')) {
    dat <- dat[!(chromosome == snakemake@config$loci[[x]]$chrom & base_pair_location %between% c(snakemake@config$loci[[x]]$start, snakemake@config$loci[[x]]$stop))]
  }
}

dat[, `:=` (Z.A = beta_a/se_a, Z.B = beta_b/se_b), env = list(beta_a = beta_a_col, beta_b = beta_b_col, se_a = se_a_col, se_b = se_b_col)]

ggsave(ggplot(dat)+
geom_point(aes(x = Z.A, y = Z.B), alpha = 0.25)+
geom_hline(yintercept = qnorm(2.5e-8), linetype = 'dashed')+
geom_hline(yintercept = -qnorm(2.5e-8), linetype = 'dashed')+
geom_vline(xintercept = qnorm(2.5e-8), linetype = 'dashed')+
geom_vline(xintercept = -qnorm(2.5e-8), linetype = 'dashed')+
xlab(snakemake@wildcards$study_a)+
ylab(snakemake@wildcards$study_b)+
xlim(c(-20, 20))+
ylim(c(-20, 20))+
coord_fixed(),
file = snakemake@output[[1]],
width = 6,
height = 6
)
