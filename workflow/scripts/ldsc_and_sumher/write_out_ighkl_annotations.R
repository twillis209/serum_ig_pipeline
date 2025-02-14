library(data.table)
setDTthreads(snakemake@threads)

snps <- fread(snakemake@input[[1]], header = F)
names(snps) <- 'ID'
snps[, c('chr', 'pos') := tstrsplit(ID, split = ':', keep = 1:2)]
snps[, `:=` (pos = as.integer(pos))]

for(x in c('igh', 'igk', 'igl')) {
  snps[chr == snakemake@config$loci[[x]]$chrom & pos %between% c(snakemake@config$loci[[x]]$start, snakemake@config$loci[[x]]$stop), set := TRUE, env = list(set = x)]
  fwrite(snps[get(x) == TRUE, .(ID)], sep = ' ', snakemake@output[[x]])
}

fwrite(unique(snps[is.na(igh)& is.na(igk) & is.na(igl), .(ID)]), sep = ' ', file = snakemake@output$base)








