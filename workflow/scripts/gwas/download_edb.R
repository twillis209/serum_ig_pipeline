library(AnnotationHub)
library(ensembldb)

ah <- AnnotationHub()

q <- query(ah, snakemake@params$search_tokens)

q <- q[q$genome == snakemake@params$genome]

if(length(q) != 1) {
  stop(sprintf("%d hits for AnnotationHub query '%s' for genome '%s'. Expected 1 hit.", length(q), paste(snakemake@params$search_tokens, collapse = ','), snakemake@params$genome))
}

edb <- q[[1]]

file.copy(edb@ensdb@dbname, snakemake@output[[1]])
## file.copy(edb@ensdb@dbname, "resources/gwas/ensembl_113_hsapiens_edb.sqlite")
