library(data.table)
setDTthreads(snakemake@threads)

iso <- snakemake@wildcards$isotype
non_ig <- snakemake@wildcards$non_ig

merged <- fread(snakemake@input$merged)

leads <- fread(snakemake@input$lead_snps)

leads <- leads[rsid_a == snakemake@wildcards$first_rsid & rsid_b == snakemake@wildcards$second_rsid, env = list(rsid_a = paste('rsid', iso, sep = '.'), rsid_b = paste('rsid', non_ig, sep = '.'))]
if(leads[, .N] == 0) {
  stop("No lead SNPs matching rsids")
}

max_bp <- leads[, max(bp_first, bp_second), env = list(bp_first = paste('base_pair_location', iso, sep = '.'), bp_second = paste('base_pair_location', non_ig, sep = '.'))]
min_bp <- leads[, min(bp_first, bp_second), env = list(bp_first = paste('base_pair_location', iso, sep = '.'), bp_second = paste('base_pair_location', non_ig, sep = '.'))]
chrom <- leads[, chrom, env = list(chrom = sprintf('chromosome.%s', iso))]

fwrite(merged[chromosome == as.character(chrom) & base_pair_location %between% c(min_bp - snakemake@params$flank, max_bp + snakemake@params$flank)], sep = '\t', file = snakemake@output[[1]])
