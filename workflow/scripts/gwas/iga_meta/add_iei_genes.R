library(data.table)

lead <- fread(snakemake@input$lead)

ieis <- fread(snakemake@input$ieis)

ieis <- na.omit(ieis, c('start_position', 'end_position', 'chromosome_name', 'hgnc_symbol'))

flank <- snakemake@params$flank

ieis[, `:=` (f_start = start_position - flank, f_end = end_position + flank, chromosome = as.character(chromosome_name))]
setkey(ieis, chromosome, f_start, f_end)

lead[, `:=` (f_start = base_pair_location, f_end = base_pair_location + 1, chromosome = as.character(chromosome))]
setkey(lead, chromosome, f_start, f_end)

merged <- foverlaps(lead, ieis, mult = 'all')[!is.na(hgnc_symbol)]

merged[base_pair_location >= start_position & base_pair_location <= end_position, distance_to_iei_gene := 0]
merged[base_pair_location <= start_position, distance_to_iei_gene := start_position - base_pair_location]
merged[base_pair_location >= end_position, distance_to_iei_gene := base_pair_location - end_position]

merged[, c('f_start', 'f_end', 'i.f_start', 'i.f_end', 'chromosome_name', 'start_position', 'end_position', 'sanitised_gene') := NULL]
setnames(merged, c('hgnc_symbol', 'disorder'), c('iei_hgnc_symbol', 'IEIs'))

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
