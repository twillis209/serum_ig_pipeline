library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col

dat <- fread(snakemake@input[[1]])

for(x in snakemake@params$loci_to_drop) {
  dat <- dat[!(chr == snakemake@config$loci[[x]]$chrom &
               pos %between% c(snakemake@config$loci[[x]]$start, snakemake@config$loci[[x]]$stop)
  ), env = list(chr = chr_col, pos = bp_col)]
}

fwrite(dat, sep = '\t', file = snakemake@output[[1]])
