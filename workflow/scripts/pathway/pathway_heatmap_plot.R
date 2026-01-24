library(data.table)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(magrittr)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Helper function to get gene symbols from Entrez IDs
get_gene_symbols <- function(entrezid, entrez_db) {
  entrez_db[entrez_db[, "ENTREZID"] %in% entrezid, "SYMBOL"]
}

# Run KEGG pathway enrichment
run_kk <- function(entrez_db) {
  kk <- enrichKEGG(
    gene = entrez_db[, "ENTREZID"],
    organism = "hsa",
    qvalueCutoff = 0.1
  )
  invisible(kk)
}

# Run Reactome pathway enrichment
run_react <- function(entrez_db) {
  rct <- enrichPathway(
    gene = entrez_db[, "ENTREZID"],
    pvalueCutoff = 0.05, 
    readable = TRUE
  )
  invisible(rct)
}

# Get Entrez database from gene symbols and Ensembl IDs
get_entrez_db <- function(xm) {
  genes <- strsplit(xm$genes, ',') %>%
    unlist() %>%
    unique()

  symbols <- genes[!grepl("^ENSG", genes)]
  ensembl_ids <- genes[grepl("^ENSG", genes)]

  res <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    columns = c("ENTREZID", "SYMBOL", "ENSEMBL"),
    keytype = "SYMBOL"
  )

  res2 <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ensembl_ids,
    columns = c("ENTREZID", "SYMBOL", "ENSEMBL"),
    keytype = "ENSEMBL"
  )

  rbind(subset(res, !is.na(ENTREZID)), subset(res2, !is.na(ENTREZID))) %>%
    unique()
}

# Function to create matrix mapping variants to pathways
getM <- function(kk, trans = TRUE) {
  dt <- as.data.table(kk)[qvalue < 0.05]
  if (trans) {
    entrez_db <- get_entrez_db(x)
    memgenes <- strsplit(dt$geneID, "/") %>%
      lapply(., get_gene_symbols, entrez_db = entrez_db)
  } else {
    memgenes <- strsplit(dt$geneID, "/")
  }
  allgenes <- unlist(memgenes) %>% unique()
  # maps variants 2 genes
  G <- matrix(0, nrow(x), length(allgenes), dimnames = list(x$rsID, allgenes))
  Gi <- strsplit(gsub(" ", "", x$genes), ",")
  for (i in 1:nrow(x)) {
    G[i, ] <- as.numeric(allgenes %in% Gi[[i]])
  }
  # maps genes 2 pathways
  P <- matrix(0, length(allgenes), length(memgenes), dimnames = list(allgenes, dt$Description))
  for (i in seq_along(memgenes)) {
    P[, i] <- as.numeric(allgenes %in% memgenes[[i]])
  }
  # combine
  M <- G %*% P  # will count through how many genes a variant is linked to a pathway
  M <- pmin(M, 1)  # binarises this
  M <- M[rowSums(M) > 0, colSums(M) > 4]
  M
}

drop_ensid_genes <- function(gene_string) {
  genes <- trimws(unlist(strsplit(gene_string, ",")))
  filtered <- genes[!grepl("^ENSG", genes)]   # use "^ENSG" to match names that start with ENSG
  result <- paste(filtered, collapse = ", ")
  result
}

# Function to plot heatmap of pathways (columns) vs variants (rows)
draw_pathway_heatmap <- function(kk, trans = TRUE, ...) {
  M <- getM(kk, trans = trans)
  on <- colnames(M)
  cann <- x[rsID %in% rownames(M), .(rsID, IgA = grepl("IgA", what), IgG = grepl("IgG", what), IgM = grepl("IgM", what))] %>%
    as.data.frame()
  cann$IgA %<>% as.numeric()
  cann$IgM %<>% as.numeric()
  cann$IgG %<>% as.numeric()
  rownames(cann) <- cann$rsID
  cann <- cann[, c("IgA", "IgG", "IgM")]

  om <- rownames(M)
  summary(nchar(om))

  m <- sprintf("%-11s\t%-22s", sub(":.*", "", om), x$genes[match(om, x$rsID)])
  m <- sapply(m, drop_ensid_genes)
  rownames(M) <- m
  rownames(cann) <- m
  ## print(dim(M))
  pheatmap(M,
    show_rownames = TRUE, 
    show_colnames = TRUE,
    color = colorRampPalette(c("grey90", "grey20"))(50),
    clustering_method = "ward.D",
    labels_col = stringr::str_wrap(colnames(M), 50, exdent = 6),
    annotation_row = cann,
    cutree_row = 7, 
    cutree_col = 5,
    legend = FALSE, 
    annotation_legend = FALSE, 
    ...
  )
}

# ============================================================================
# Main analysis workflow
# ============================================================================

# Read input data from snakemake
x <- fread(snakemake@input[[1]])

# Get Entrez database
entrez <- get_entrez_db(x)

# Run KEGG enrichment
kk <- run_kk(entrez)

# Run Reactome enrichment
rct <- run_react(entrez)

# Create combined Reactome + KEGG heatmap
memgenes <- strsplit(kk$geneID, "/") %>%
  lapply(., get_gene_symbols, entrez_db = entrez) %>%
  sapply(., paste, collapse = "/")
kk2 <- as.data.table(kk)
kk2$geneID <- memgenes
m <- rbind(kk2, as.data.table(rct), fill = TRUE)

height <- 14
width <- 12
cellwidth <- 13
cellheight <- 10
fontsize_row <- 10
fontsize_col <- 8

for (x in c(snakemake@output[snakemake@output != ""])) {
  draw_pathway_heatmap(m, FALSE,
    filename = x,
    width = width,
    height = height,
    cellwidth = cellwidth,
    cellheight = cellheight,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col
  )
}
