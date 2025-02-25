library(data.table)

lead <- fread(snakemake@input$lead)

ieis <- fread(snakemake@input$ieis, select = c('symbol', 'chromosome', 'start', 'end'))

flank <- 2e5

ieis[, `:=` (f_start = start - flank, f_end = end + flank, chromosome = as.character(chromosome))]
setkey(ieis, chromosome, f_start, f_end)

lead[, `:=` (f_start = base_pair_location, f_end = base_pair_location + 1, chromosome = as.character(chromosome))]
setkey(lead, chromosome, f_start, f_end)

merged <- foverlaps(lead, ieis, mult = 'all')[!is.na(symbol)]

merged[base_pair_location >= start & base_pair_location <= end, distance_to_iei_gene := 0]
merged[base_pair_location <= start, distance_to_iei_gene := start - base_pair_location]
merged[base_pair_location >= end, distance_to_iei_gene := base_pair_location - end]

merged[, c('f_start', 'f_end', 'i.f_start', 'i.f_end') := NULL]
setnames(merged, 'symbol', 'iei_gene')

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
