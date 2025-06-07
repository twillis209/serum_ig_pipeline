library(data.table)
setDTthreads(snakemake@threads)
library(coloc)

beta_a <- sprintf('beta.%s', snakemake@wildcards$first_isotype)
beta_b <- sprintf('beta.%s', snakemake@wildcards$second_isotype)
se_a <- sprintf('standard_error.%s', snakemake@wildcards$first_isotype)
se_b <- sprintf('standard_error.%s', snakemake@wildcards$second_isotype)
z_a <- sprintf('z.%s', snakemake@wildcards$first_isotype)
z_b <- sprintf('z.%s', snakemake@wildcards$second_isotype)
n_a <- sprintf('sample_size.%s', snakemake@wildcards$first_isotype)
n_b <- sprintf('sample_size.%s', snakemake@wildcards$second_isotype)

dat <- fread(snakemake@input[[1]])

dat <- unique(dat, by = 'rsid')

dat <- na.omit(dat, c(beta_a, beta_b, se_a, se_b))

if(snakemake@wildcards$trim == 'trimmed') {
  dat[, n_frac_a := as.integer(n_a)/as.integer(snakemake@params$first_isotype_max_n), env = list(n_a = n_a)]
  dat[, n_frac_b := as.integer(n_b)/as.integer(snakemake@params$second_isotype_max_n), env = list(n_b = n_b)]

  dat <- dat[abs(n_frac_a - n_frac_b) < 0.1]
}

min_p_first_isotype <- dat[, min(p), env = list(p = sprintf("p_value.%s", snakemake@wildcards$first_isotype))]
min_p_second_isotype <- dat[, min(p), env = list(p = sprintf("p_value.%s", snakemake@wildcards$second_isotype))]

first_dataset <- list(snp = dat[, rsid],
                      beta = dat[, beta, env = list(beta = beta_a)],
                      varbeta = dat[, se^2, env = list(se = se_a)],
                      type = 'quant',
                      N = dat[, n, env = list(n = n_a)],
                      sdY = 1
                      )

second_dataset <- list(snp = dat[, rsid],
                      beta = dat[, beta, env = list(beta = beta_b)],
                      varbeta = dat[, se^2, env = list(se = se_b)],
                      type = 'quant',
                      N = dat[, n, env = list(n = n_b)],
                      sdY = 1
                      )

res <- coloc.abf(first_dataset, second_dataset)

dat[, `:=` (z_a = beta_a/se_a, z_b = beta_b/se_b), env = list(z_a = z_a, beta_a = beta_a, se_a = se_a, z_b = z_b, beta_b = beta_b, se_b = se_b)]

z_cor <- dat[, cor(z_a, z_b, use = "pairwise.complete.obs", method = "pearson"),
             env = list(z_a = z_a, z_b = z_b)]

saveRDS(res, file = snakemake@output[['rds']])

# Ratio of effect estimates at lead SNPs
first_iso_lead_snp_effect_ratio <- dat[rsid == snakemake@wildcards$first_rsid, beta_a/beta_b, env = list(beta_a = beta_a, beta_b = beta_b)]
second_iso_lead_snp_effect_ratio <- dat[rsid == snakemake@wildcards$second_rsid, beta_a/beta_b, env = list(beta_a = beta_a, beta_b = beta_b)]

res_dat <- data.table(t(res$summary))

res_dat[, `:=` (first_trait = snakemake@wildcards$first_isotype, second_trait = snakemake@wildcards$second_isotype, trimmed = snakemake@wildcards$trim == 'trimmed', first_snp = snakemake@wildcards$first_rsid, second_snp = snakemake@wildcards$second_rsid, min_p.first = min_p_first_isotype, min_p.second = min_p_second_isotype, pearson.cor = z_cor, first_iso_lead_snp_effect_ratio = first_iso_lead_snp_effect_ratio, second_iso_lead_snp_effect_ratio = second_iso_lead_snp_effect_ratio)]

fwrite(res_dat, file = snakemake@output$tsv, sep = '\t')
