library(data.table)
setDTthreads(snakemake@threads)
library(coloc)

#save.image('coloc.RData')

beta_a <- sprintf('beta.%s', snakemake@wildcards$first_isotype)
beta_b <- sprintf('beta.%s', snakemake@wildcards$second_isotype)
se_a <- sprintf('standard_error.%s', snakemake@wildcards$first_isotype)
se_b <- sprintf('standard_error.%s', snakemake@wildcards$second_isotype)
n_a <- sprintf('sample_size.%s', snakemake@wildcards$first_isotype)
n_b <- sprintf('sample_size.%s', snakemake@wildcards$second_isotype)

dat <- fread(snakemake@input[[1]])

dat <- unique(dat, by = 'rsid')

dat <- na.omit(dat, c(beta_a, beta_b, se_a, se_b))

dat[, n_frac_a := as.integer(n_a)/as.integer(snakemake@params$first_isotype_max_n), env = list(n_a = n_a)]
dat[, n_frac_b := as.integer(n_b)/as.integer(snakemake@params$second_isotype_max_n), env = list(n_b = n_b)]

dat <- dat[abs(n_frac_a - n_frac_b) < 0.1]

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

saveRDS(coloc.abf(first_dataset, second_dataset), file = snakemake@output[[1]])
