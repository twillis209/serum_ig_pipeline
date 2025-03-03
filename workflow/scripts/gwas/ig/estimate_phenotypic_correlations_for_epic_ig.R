library(data.table)
library(boot)

dat <- fread(snakemake@input[[1]])

dat[, IGG_CONC := IGG1_CONC + IGG2_CONC + IGG3_CONC + IGG4_CONC]

cor_fun <- function(dtbl, indices, iso_a, iso_b, transform = identity) dtbl[indices, cor(transform(iso_a), transform(iso_b), use = 'pairwise.complete.obs'), env = list(iso_a = iso_a, iso_b = iso_b)]

set.seed(snakemake@params$seed)

combs <- combn(c('IGA_CONC', 'IGG_CONC', 'IGM_CONC'), 2)

est <- numeric()
est_lower_bound <- numeric()
est_upper_bound <- numeric()

 for(i in 1:ncol(combs)) {
  raw_res <- boot(dat, statistic = cor_fun, R = snakemake@params$reps, iso_a = combs[1,i], iso_b = combs[2,i])

  raw_ci <- boot.ci(raw_res, type = 'basic')

  est <- c(est, raw_ci$t0)
  est_lower_bound <- c(est_lower_bound, raw_ci$basic[,4])
  est_upper_bound <- c(est_upper_bound, raw_ci$basic[,5])

  log_res <- boot(dat, statistic = cor_fun, R = snakemake@params$reps, iso_a = combs[1,i], iso_b = combs[2,i], transform = log)

  log_ci <- boot.ci(log_res, type = 'basic')

  est <- c(est, log_ci$t0)
  est_lower_bound <- c(est_lower_bound, log_ci$basic[,4])
  est_upper_bound <- c(est_upper_bound, log_ci$basic[,5])
}

res_dat <- data.table(scale = rep(c('raw', 'log'), 3),
           isotype_a = rep(combs[1,], each = 2),
           isotype_b = rep(combs[2,], each = 2),
           est,
           est_lower_bound,
           est_upper_bound
           )
setkey(res_dat, isotype_a)

name_dat <- data.table(raw = c('IGA_CONC', 'IGG_CONC', 'IGM_CONC'), pretty = c('IgA', 'IgG', 'IgM'))
setkey(name_dat, raw)

res_dat[name_dat, isotype_a := pretty]
setkey(res_dat, isotype_b)
res_dat[name_dat, isotype_b := pretty]

fwrite(res_dat[, .(`First isotype` = isotype_a,
                   `Second isotype` = isotype_b,
                   Scale = scale,
                   Estimate = est,
                   `95% CI lower bound` = est_lower_bound,
                   `95% CI upper bound` = est_upper_bound
                   )],
       file = snakemake@output[[1]], sep = '\t')
