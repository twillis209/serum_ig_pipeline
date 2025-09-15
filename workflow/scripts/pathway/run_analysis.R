library(data.table)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(ggtext)
library(stringr)
library(patchwork)
library(cowplot)
library(dplyr)
library(magrittr)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)

get_gene_symbols <- function(entrezid, entrez_db) {
  entrez_db[entrez_db[, "ENTREZID"] %in% entrezid, "SYMBOL"]
}

run_kk <- function(entrez_db) {
  kk <- enrichKEGG(
    gene = entrez_db[, "ENTREZID"],
    organism = "hsa",
    qvalueCutoff = 0.1
  )
  invisible(kk)
}
run_react <- function(entrez_db) {
  rct <- enrichPathway(
    gene = entrez_db[, "ENTREZID"],
    pvalueCutoff = 0.05, readable = TRUE
  )
  invisible(rct)
}

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


  ## entrez_db <- bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

  ## entrez_db

  rbind(subset(res, !is.na(ENTREZID)), subset(res2, !is.na(ENTREZID))) %>%
    unique()
}

x <- fread(snakemake@input[[1]])

kk <- run_kk(get_entrez_db(x))

rct <- run_react(get_entrez_db(x))

## make data structure and group variants by pathway analysis
getM <- function(kk, trans = TRUE) {
  dt <- as.data.table(kk)[qvalue < 0.05]
  if (trans) {
    entrez_db <- get_entrez_db(x)
    memgenes <- strsplit(dt$geneID, "/") %>%
      lapply(., get_gene_symbols)
  } else {
    memgenes <- strsplit(dt$geneID, "/")
  }
  allgenes <- unlist(memgenes) %>% unique()
  ## maps variants 2 genes
  G <- matrix(0, nrow(x), length(allgenes), dimnames = list(x$rsID, allgenes))
  Gi <- strsplit(gsub(" ", "", x$genes), ",")
  for (i in 1:nrow(x)) {
    G[i, ] <- as.numeric(allgenes %in% Gi[[i]])
  }
  ## maps genes 2 pathways
  P <- matrix(0, length(allgenes), length(memgenes), dimnames = list(allgenes, dt$Description))
  for (i in seq_along(memgenes)) {
    P[, i] <- as.numeric(allgenes %in% memgenes[[i]])
  }
  ## combine
  M <- G %*% P # will count through how many genes a variant is linked to a pathway
  M <- pmin(M, 1) # binarises this
  M <- M[rowSums(M) > 0, colSums(M) > 4]
  M
}

# Draws heatmaps of the pathway enrichment results
draw_pathway_heatmap <- function(kk, trans = TRUE, ...) {
  M <- getM(kk, trans = trans)
  # M=t(M)
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
  rownames(M) <- m
  rownames(cann) <- m
  print(dim(M))
  pheatmap(M,
    show_rownames = TRUE, show_colnames = TRUE,
    color = colorRampPalette(c("grey90", "grey20"))(50),
    clustering_method = "ward.D",
    ## labels_col=substring(colnames(M),1,nchar(colnames(M))),
    labels_col = stringr::str_wrap(colnames(M), 50, exdent = 6),
    annotation_row = cann,
    # annotation_row=rann,
    cutree_row = 7, cutree_col = 5,
    legend = FALSE, annotation_legend = FALSE, ...
  )
}

entrez <- get_entrez_db(x)

memgenes <- strsplit(kk$geneID, "/") %>%
  lapply(., get_gene_symbols, entrez_db = entrez) %>%
  sapply(., paste, collapse = "/")
kk2 <- as.data.table(kk)
kk2$geneID <- memgenes
m <- rbind(kk2, as.data.table(rct), fill = TRUE)

draw_pathway_heatmap(m, FALSE, main = "Reactome + KEGG", filename = snakemake@output$combined_pathways_heatmap, width = 12, height = 10, fontsize = 8)

kk_for_plot <- as.data.table(kk)[qvalue < 0.05]
kk_for_plot <- kk_for_plot[order(-qvalue)]
kk_for_plot[, Description := factor(Description, levels = Description)]
cols <- enrichplot:::get_enrichplot_color(2)
kk_pl <- ggplot(kk_for_plot, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment), fill = -log10(qvalue), size = Count)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = cols[2], high = cols[1]) +
  scale_x_continuous(limits = c(1, max(kk_for_plot$FoldEnrichment) + .5)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20)) +
  theme_cowplot() +
  labs(y = "", x = "Fold enrichment") +
  guides(fill = guide_colorbar(title = "-log10(q-value)"))+
  background_grid(major = "y") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 14))

rct_for_plot <- as.data.table(rct)[qvalue < 0.05]
rct_for_plot <- rct_for_plot[order(-qvalue)]
rct_for_plot[, Description := factor(Description, levels = Description)]
rct_for_plot[, label_len := nchar(as.character(Description))]
rct_for_plot[, label := gsub("\n", "<br>", sprintf("<span style='font-size:14pt;'>%s</span>", str_wrap(Description, width = 50)))]

cols <- enrichplot:::get_enrichplot_color(2)
rct_pl <- ggplot(rct_for_plot, aes(x = FoldEnrichment, y = reorder(label, FoldEnrichment), fill = -log10(qvalue), size = Count)) +
  geom_point(shape = 21) +
  scale_fill_gradient(low = cols[2], high = cols[1]) +
  scale_x_continuous(limits = c(1, max(rct_for_plot$FoldEnrichment) + .5)) +
  theme_cowplot() +
  labs(y = "", x = "Fold enrichment") +
  background_grid(major = "y") +
  theme(axis.text.y = element_markdown())+
  theme(legend.position = "none")

enrichment_pls <- (kk_pl | rct_pl) + plot_annotation(tag_levels = "A")

ggsave(enrichment_pls, filename = snakemake@output$combined_pathways_enrichment_plot, width = 12, height = 12)
