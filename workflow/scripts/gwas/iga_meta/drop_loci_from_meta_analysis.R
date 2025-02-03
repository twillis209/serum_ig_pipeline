library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col

dat <- fread(snakemake@input[[1]])

# IGH
dat <- dat[!(chr == snakemake@config$loci$igh$chrom &
             pos %between% c(snakemake@config$loci$igh$start, snakemake@config$loci$igh$stop)
            ), env = list(chr = chr_col, pos = bp_col)]

fwrite(dat, sep = '\t', file = snakemake@output[[1]])
