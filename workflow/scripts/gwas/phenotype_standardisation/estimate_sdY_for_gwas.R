library(data.table)
setDTthreads(snakemake@threads)
library(tomics)

se_col <- snakemake@config$se_col

dat <- fread(snakemake@input[[1]])

dat <- dat[!is.na(ALT_FREQS)]

sdY.estimates <- numeric(snakemake@params$reps)

set.seed(snakemake@params$seed)

for(i in 1:snakemake@params$reps) {
  sample_indices <- sample(1:nrow(dat), size = snakemake@params$no_of_sample_variants)
  sdY.estimates[i] <- sdY.est(vbeta = dat[sample_indices, se_col^2, env = list(se_col = se_col)], maf = dat[sample_indices, ALT_FREQS], n = snakemake@params$N)
  dat <- dat[-sample_indices]
}

fwrite(data.table(sdY.est = sdY.estimates), file = snakemake@output[[1]], sep = '\t')
