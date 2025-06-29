library(data.table)

dat <- fread(snakemake@input, sep = '\t', header = T)

