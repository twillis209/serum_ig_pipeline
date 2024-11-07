library(data.table)

iga <- fread(snakemake@input$iga)
iei <- fread(snakemake@input$iei)

iga[, `:=` (bp38.start = bp38, bp38.stop = bp38 + 1)]
setkey(iga, chr38, bp38.start, bp38.stop)

setnames(iei, c('start', 'end'), c('start.bp', 'end.bp'))
iei[, `:=` (start.50kb = start.bp - 5e4, end.50kb = end.bp + 5e4)]
setkey(iei, chromosome, start.50kb, end.50kb)

overlaps <- foverlaps(iga, iei, mult = 'all', nomatch = NA)[!is.na(start.50kb)]

fwrite(overlaps[, .(gene_name, genes, chr19, start.bp, end.bp, bp38.start, rsID)], sep = '\t', file = snakemake@output[[1]])
