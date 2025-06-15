library(data.table)
setDTthreads(snakemake@threads)
library(coloc)

beta_a <- sprintf('beta.%s', snakemake@wildcards$isotype)
beta_b <- sprintf('beta.%s', snakemake@wildcards$non_ig)
z_a <- sprintf('z.%s', snakemake@wildcards$isotype)
z_b <- sprintf('z.%s', snakemake@wildcards$non_ig)
se_a <- sprintf('standard_error.%s', snakemake@wildcards$isotype)
se_b <- sprintf('standard_error.%s', snakemake@wildcards$non_ig)
n_a <- sprintf('sample_size.%s', snakemake@wildcards$isotype)
n_b <- sprintf('sample_size.%s', snakemake@wildcards$non_ig)

dat <- fread(snakemake@input[[1]])

dat <- unique(dat, by = 'rsid')

dat <- na.omit(dat, c(beta_a, beta_b, se_a, se_b))

min_p_ig <- dat[, min(p), env = list(p = sprintf("p_value.%s", snakemake@wildcards$isotype))]
min_p_non_ig <- dat[, min(p), env = list(p = sprintf("p_value.%s", snakemake@wildcards$non_ig))]
min_p_rsid_non_ig <- dat[which.min(p), rsid, env = list(p = sprintf("p_value.%s", snakemake@wildcards$non_ig))]

ig_dataset <- list(snp = dat[, rsid],
                      beta = dat[, beta, env = list(beta = beta_a)],
                      varbeta = dat[, se^2, env = list(se = se_a)],
                      type = 'quant',
                      N = snakemake@config$gwas_datasets[[snakemake@wildcards$isotype]]$samples,
                      sdY = 1
                   )

if(snakemake@wildcards$non_ig == 'lymphocyte-counts') {
  non_ig_dataset <- list(snp = dat[, rsid],
                        beta = dat[, beta, env = list(beta = beta_b)],
                        varbeta = dat[, se^2, env = list(se = se_b)],
                        type = 'quant',
                        N = snakemake@config$gwas_datasets[[snakemake@wildcards$non_ig]]$samples,
                        sdY = 1
                        )
} else {
  non_ig_dataset <- list(snp = dat[, rsid],
                        beta = dat[, beta, env = list(beta = beta_b)],
                        varbeta = dat[, se^2, env = list(se = se_b)],
                        type = 'cc',
                        N = snakemake@config$gwas_datasets[[snakemake@wildcards$non_ig]]$samples,
                        s = snakemake@config$gwas_datasets[[snakemake@wildcards$non_ig]]$cases/snakemake@config$gwas_datasets[[snakemake@wildcards$non_ig]]$samples,
                        sdY = 1
                        )
}

res <- coloc.abf(ig_dataset, non_ig_dataset)

dat[, `:=` (z_a = beta_a/se_a, z_b = beta_b/se_b), env = list(z_a = z_a, beta_a = beta_a, se_a = se_a, z_b = z_b, beta_b = beta_b, se_b = se_b)]

z_cor <- dat[, cor(z_a, z_b, use = "pairwise.complete.obs", method = "pearson"),
             env = list(z_a = z_a, z_b = z_b)]

saveRDS(res, file = snakemake@output[['rds']])

# Ratio of effect estimates at lead SNPs
ig_lead_snp_effect_ratio <- dat[rsid == snakemake@wildcards$isotype_rsid, beta_a/beta_b, env = list(beta_a = beta_a, beta_b = beta_b)]
non_ig_lead_snp_effect_ratio <- dat[rsid == min_p_rsid_non_ig, beta_a/beta_b, env = list(beta_a = beta_a, beta_b = beta_b)]

res_dat <- data.table(t(res$summary))

res_dat[,
 `:=` (first_trait = snakemake@wildcards$isotype,
 second_trait = snakemake@wildcards$non_ig,
 ig_snp = snakemake@wildcards$isotype_rsid,
 non_ig_snp = min_p_rsid_non_ig,
 chromosome = dat[rsid == snakemake@wildcards$isotype_rsid, chromosome],
 ig_snp_pos = dat[rsid == snakemake@wildcards$isotype_rsid, base_pair_location],
 non_ig_snp_pos = dat[rsid == min_p_rsid_non_ig, base_pair_location],
 min_p.first = min_p_ig,
 min_p.second = min_p_non_ig,
 pearson.cor = z_cor,
 ig_snp_effect_ratio = ig_lead_snp_effect_ratio,
 non_ig_snp_effect_ratio = non_ig_lead_snp_effect_ratio)]

fwrite(res_dat, file = snakemake@output$tsv, sep = '\t')
