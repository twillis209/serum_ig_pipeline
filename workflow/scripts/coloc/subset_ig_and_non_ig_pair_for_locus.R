library(data.table)
setDTthreads(snakemake@threads)

iso <- snakemake@wildcards$isotype

merged <- fread(snakemake@input$merged)

leads <- fread(snakemake@input$lead_snps)

leads <- leads[rsid == snakemake@wildcards$isotype_rsid]
if(leads[, .N] == 0) {
  stop("No lead SNPs matching rsids")
}

chrom <- leads[, chromosome]
min_bp <- max(leads[, base_pair_location] - snakemake@params$flank, 0)
max_bp <- leads[, base_pair_location] + snakemake@params$flank

names(merged) <- gsub('-meta', '', names(merged))

fwrite(merged[chromosome == as.character(chrom) & base_pair_location %between% c(min_bp, max_bp)], sep = '\t', file = snakemake@output[[1]])
