library(data.table)
setDTthreads(snakemake@threads)
#save.image('sdY.RData')

dat <- fread(snakemake@input[[1]])

beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col

dat[, beta_col := beta_col/snakemake@params$sdY_estimate, env = list(beta_col = beta_col)]
dat[, se_col := se_col/snakemake@params$sdY_estimate, env = list(se_col = se_col)]

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
