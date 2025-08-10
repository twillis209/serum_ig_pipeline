library(GenomicRanges)
library(ensembldb)

edb <- EnsDb(snakemake@input$edb)

daf <- read.table(snakemake@input$lead, header = TRUE, sep = '\t')

daf$chromosome <- as.character(daf$chromosome)

daf[daf$chromosome == '23', 'chromosome'] <- 'X'

snps_gr <- GRanges(
  seqnames = daf$chromosome,
  ranges = IRanges(start = daf$base_pair_location, width = 1),
  rsid = daf$rsid
)

genes_gr <- genes(edb, filter = GeneBiotypeFilter("protein_coding"))

suppressWarnings(nearest_hits <- distanceToNearest(snps_gr, genes_gr))

nearest_genes <- data.frame(
  rsid = mcols(snps_gr[queryHits(nearest_hits)])$rsid,
  nearest_gene_id        = names(genes_gr[subjectHits(nearest_hits)]),
  nearest_gene_name      = genes_gr[subjectHits(nearest_hits)]$gene_name,
  distance_bp    = mcols(nearest_hits)$distance
)

merged <- merge(daf, nearest_genes, by = 'rsid', all.x = TRUE)

merged[merged$chromosome == 'X', 'chromosome'] <- '23'

write.table(merged, file = snakemake@output[[1]], sep = '\t', row.names = FALSE, quote = FALSE)
