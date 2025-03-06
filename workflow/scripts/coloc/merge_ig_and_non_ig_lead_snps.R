library(data.table)

cols <- c('chromosome', 'rsid', 'base_pair_location')

flank <- snakemake@params$window/2

iso <- snakemake@wildcards$isotype
non_ig <- snakemake@wildcards$non_ig

ig <- fread(snakemake@input$ig, select = cols)
setnames(ig, cols, paste(cols, iso, sep = '.'))
ig[, `:=` (start = bp - flank, end = bp + flank), env = list(start = sprintf("start.%s", iso),
                                                            bp = sprintf("base_pair_location.%s", iso),
                                                            end = sprintf("end.%s", iso)
                                                            )
   ]
setkeyv(ig, c(sprintf("chromosome.%s", iso), sprintf("start.%s", iso), sprintf("end.%s", iso)))

non_ig_dat <- fread(snakemake@input$non_ig, select = cols)
setnames(non_ig_dat, cols, paste(cols, non_ig, sep = '.'))
non_ig_dat[, `:=` (start = bp - flank, end = bp + flank), env = list(start = sprintf("start.%s", non_ig),
                                                            bp = sprintf("base_pair_location.%s", non_ig),
                                                            end = sprintf("end.%s", non_ig)
                                                            )
   ]
setkeyv(non_ig_dat, c(sprintf("chromosome.%s", non_ig), sprintf("start.%s", non_ig), sprintf("end.%s", non_ig)))

ig_non_ig <- foverlaps(ig, non_ig_dat, mult = 'all')[!is.na(rsid), env = list(rsid = sprintf("rsid.%s", iso))]
ig_non_ig[, distance := abs(bp1 - bp2), env = list(bp1 = sprintf("base_pair_location.%s", iso), bp2 = sprintf("base_pair_location.%s", non_ig))]
fwrite(ig_non_ig, file = snakemake@output[[1]], sep = '\t')
