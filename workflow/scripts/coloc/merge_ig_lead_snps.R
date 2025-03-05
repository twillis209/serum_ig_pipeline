library(data.table)

cols <- c('chromosome', 'rsid', 'base_pair_location')

flank <- snakemake@params$window/2

iga <- fread(snakemake@input$iga, select = cols)
setnames(iga, cols, paste(cols, 'iga', sep = '.'))
iga[, `:=` (start.iga = base_pair_location.iga - flank, end.iga = base_pair_location.iga + flank)]
setkey(iga, chromosome.iga, start.iga, end.iga)

igg <- fread(snakemake@input$igg, select = cols)
setnames(igg, cols, paste(cols, 'igg', sep = '.'))
igg[, `:=` (start.igg = base_pair_location.igg - flank, end.igg = base_pair_location.igg + flank)]
setkey(igg, chromosome.igg, start.igg, end.igg)

igm <- fread(snakemake@input$igm, select = cols)
setnames(igm, cols, paste(cols, 'igm', sep = '.'))
igm[, `:=` (start.igm = base_pair_location.igm - flank, end.igm = base_pair_location.igm + flank)]
setkey(igm, chromosome.igm, start.igm, end.igm)

iga_igg <- foverlaps(iga, igg, mult = 'all')[!is.na(rsid.igg)]
iga_igg[, distance := abs(base_pair_location.iga - base_pair_location.igg)]
fwrite(iga_igg, file = snakemake@output$iga_igg, sep = '\t')

iga_igm <- foverlaps(iga, igm, mult = 'all')[!is.na(rsid.igm)]
iga_igm[, distance := abs(base_pair_location.iga - base_pair_location.igm)]
fwrite(iga_igm, file = snakemake@output$iga_igm, sep = '\t')

igg_igm <- foverlaps(igg, igm, mult = 'all')[!is.na(rsid.igm)]
igg_igm[, distance := abs(base_pair_location.igg - base_pair_location.igm)]
fwrite(igg_igm, file = snakemake@output$igg_igm, sep = '\t')
