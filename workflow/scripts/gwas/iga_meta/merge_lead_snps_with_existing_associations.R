library(data.table)

existing <- fread(snakemake@input$existing, select = c('chromosome', 'base_pair_location', 'rsid', 'Genes', 'dataset'))

meta <- fread(snakemake@input$meta, select = c('chromosome', 'base_pair_location', 'rsid', 'genes'))

bp_flank <- snakemake@params$window / 2

meta[, `:=` (left = base_pair_location - bp_flank,
             right = base_pair_location + bp_flank)]
setkey(meta, chromosome, left, right)

existing[, `:=` (left = base_pair_location - bp_flank,
                 right = base_pair_location + bp_flank)]
setkey(existing, chromosome, left, right)

overlaps <- foverlaps(meta, existing)

fwrite(overlaps[is.na(left), .(chromosome, base_pair_location = i.base_pair_location, genes)], sep = '\t', file = snakemake@output[[1]])
