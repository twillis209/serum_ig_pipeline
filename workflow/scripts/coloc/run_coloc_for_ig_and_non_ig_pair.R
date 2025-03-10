library(data.table)
setDTthreads(snakemake@threads)
library(coloc)

beta_a <- sprintf('beta.%s', snakemake@wildcards$isotype)
beta_b <- sprintf('beta.%s', snakemake@wildcards$non_ig)
se_a <- sprintf('standard_error.%s', snakemake@wildcards$isotype)
se_b <- sprintf('standard_error.%s', snakemake@wildcards$non_ig)
n_a <- sprintf('sample_size.%s', snakemake@wildcards$isotype)
n_b <- sprintf('sample_size.%s', snakemake@wildcards$non_ig)

dat <- fread(snakemake@input[[1]])

dat <- unique(dat, by = 'rsid')

dat <- na.omit(dat, c(beta_a, beta_b, se_a, se_b))

#dat[, n_frac_a := as.integer(n_a)/as.integer(snakemake@params$first_isotype_max_n), env = list(n_a = n_a)]
#dat[, n_frac_b := as.integer(n_b)/as.integer(snakemake@params$second_isotype_max_n), env = list(n_b = n_b)]

#dat <- dat[abs(n_frac_a - n_frac_b) < 0.1]

min_p_ig <- dat[, min(p), env = list(p = sprintf("p_value.%s", snakemake@wildcards$isotype))]
min_p_non_ig <- dat[, min(p), env = list(p = sprintf("p_value.%s", snakemake@wildcards$non_ig))]

ig_dataset <- list(snp = dat[, rsid],
                      beta = dat[, beta, env = list(beta = beta_a)],
                      varbeta = dat[, se^2, env = list(se = se_a)],
                      type = 'quant',
                      N = snakemake@config$gwas_datasets[[snakemake@wildcards$isotype]]$samples,
                      sdY = 1
                   )

if(snakemake@wildcards$non_ig == 'lymphocyte-counts') {
  non_ig_dataset <- list(snp = dat[, rsid],
                        beta = dat[, beta, env = list(beta = beta_a)],
                        varbeta = dat[, se^2, env = list(se = se_a)],
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

saveRDS(res, file = snakemake@output[['rds']])

res_dat <- data.table(t(res$summary))

res_dat[, `:=` (first_trait = snakemake@wildcards$isotype, second_trait = snakemake@wildcards$non_ig, first_snp = snakemake@wildcards$isotype_rsid, second_snp = snakemake@wildcards$non_ig_rsid, min_p.first = min_p_ig, min_p.second = min_p_non_ig)]

fwrite(res_dat, file = snakemake@output$tsv, sep = '\t')
