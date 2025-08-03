library(tomics)
library(data.table)
setDTthreads(snakemake@threads)

sumstats <- fread(snakemake@input$sumstats)
maf <- fread(snakemake@input$maf)
prune <- fread(snakemake@input$prune)

# TODO will need to split prune variant IDs so can merge
# TODO should have code to do this

# TODO old code
gwas[, SNPID := paste(CHR38, BP38, REF, ALT, sep = ':')]

merged <- merge(gwas, maf[, .(ID, ALT_FREQS)], by.x = 'SNPID', by.y = 'ID')

merged[, MAF := ifelse(ALT_FREQS > 0.5, 1-ALT_FREQS, ALT_FREQS)]

merged[, ALT_FREQS := NULL]

setkey(merged, CHR38, BP38, REF, ALT)
setkey(prune, CHR38, BP38, REF, ALT)

merged[prune, in_pruned := T]

fwrite(merged, sep = '\t', file = snakemake@output[[1]])
