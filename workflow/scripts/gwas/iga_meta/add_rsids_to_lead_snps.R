library(data.table)
setDTthreads(snakemake@threads)

lead <- fread(snakemake@input$lead)

rsids <- fread(snakemake@input$rsids)

lead[, rsid := NULL]

merged <- merge(lead, rsids, all.x = T, by = c(snakemake@config$chr_col, snakemake@config$bp_col, snakemake@config$ref_col, snakemake@config$alt_col))

fwrite(merged, sep = '\t', file = snakemake@output[[1]])
