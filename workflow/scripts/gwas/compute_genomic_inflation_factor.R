library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@params[['chr_col']]
bp_col <- snakemake@params[['bp_col']]
p_col <- snakemake@params[['p_col']]

dat <- fread(snakemake@input[[1]], sep = '\t', select = c(chr_col, bp_col, p_col))

dat[, (p_col) := as.numeric(get(p_col))]

if(!is.null(snakemake@wildcards[['variant_set']])) {
  if(snakemake@wildcards[['variant_set']] == 'sans_mhc') {
    dat <- dat[!(get(chr_col) == 6 & get(bp_col) %between% c(24e6, 45e6))]
  }
}

gif_dat <- data.table()

for(x in snakemake@params[['percentiles']]) {
  lambda <- quantile(qchisq(dat[[snakemake@params[['p_col']]]], lower.tail = F, df = 1), probs = x/100, na.rm = T)/qchisq(x/100, df = 1)

  gif_dat[, paste('lambda_0', x, sep = '_') := lambda]
}

if(!is.null(snakemake@params$controls) & !is.null(snakemake@params$cases)) {
  N0 <- as.integer(snakemake@params$controls)
  N1 <- as.integer(snakemake@params$cases)

  if(N0 == 0 & N1 > 0) {
    gif_dat[, lambda_1000 := 1  + (gif_dat[, 'lambda_0_50']-1) * (N1^-1) / (1/1000)]
  } else if(N0 > 0 & N1 == 0) {
    gif_dat[, lambda_1000 := 1 + (gif_dat[, 'lambda_0_50']-1) * (N0^-1) / (1/1000)]
  } else {
    gif_dat[, lambda_1000 := 1 + (gif_dat[, 'lambda_0_50']-1) * (N0^-1 + N1^-1) / (1/1000 + 1/1000)]
  }
}

fwrite(gif_dat, sep = '\t', file = snakemake@output[[1]])
