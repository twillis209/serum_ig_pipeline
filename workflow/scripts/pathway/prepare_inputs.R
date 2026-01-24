library(data.table)
library(dplyr)
library(magrittr)

read_and_process_ig <- function(nm) {
  x <- fread(nm)
  setnames(x, make.names(names(x)))
  setnames(x, "Gene.s.", "genes")
  x[, genes := gsub(" ", "", genes)]
  x
}

process_ig_and_coloc_results <- function(iga, igg, igm, ig_coloc) {
  iga[, what := "IgA"]
  igg[, what := "IgG"]
  igm[, what := "IgM"]
  x <- rbind(iga, igg, igm)
  x[, vid := paste0(Chromosome, "_", Position, "_", Effect.allele, "_", Other.allele)]
  x[, avid := paste0(Chromosome, "_", Position, "_", Other.allele, "_", Effect.allele)]

  ## merging colocalised Ig snps & their annotations
  setnames(ig_coloc, make.names(names(ig_coloc)))
  ig_coloc <- merge(ig_coloc, unique(x[, .(rsID, vid1 = vid)]), by.x = "First.isotype.s.lead.SNP", by.y = "rsID")
  ig_coloc <- merge(ig_coloc, unique(x[, .(rsID, vid2 = vid)]), by.x = "Second.isotype.s.lead.SNP", by.y = "rsID")
  head(ig_coloc, 2)

  ## ## decision threshold
  ## thr <- seq(0, 1, by = 0.01)
  ## fdr <- sapply(thr, function(a) with(ig_coloc[PP.H4.abf > a & Filtered == FALSE], mean(1 - PP.H4.abf)))
  ## ## plot(thr,fdr)
  ## which(fdr < 0.05) # 81 up
  ## thr[which(fdr < 0.05)] # 0.80 up
  ## ## hist(ig_coloc[Filtered==FALSE]$PP.H4.abf)
  ## ## use 0.8 to be safe
  setnames(ig_coloc, "First.isotype.s.lead.SNP", "rsid1")
  setnames(ig_coloc, "Second.isotype.s.lead.SNP", "rsid2")
  ig_coloc <- ig_coloc[PP.H4.abf >= 0.8 & Filtered == FALSE, .(vid1, vid2, rsid1, rsid2, First.isotype, Second.isotype, Pearson.correlation, PP.H4.abf)]
  ## ig_coloc # 5, all H4 > 0.93, all IgA/IgM

  ## merge rows of x
  xm <- copy(x)
  for (i in 1:nrow(ig_coloc)) {
    message(i)
    snps <- c(ig_coloc$vid1[i], ig_coloc$vid2[i])
    w <- which(x$vid %in% snps)
    if (length(w) != 2) { # unanticipated
      stop()
    }
    cat("merging ", snps[1], " ", snps[2], "\n")
    xm$vid[w[1]] <- paste0(unique(snps), collapse = ",")
    xm$what[w[1]] <- paste0(unique(xm$what[w]), collapse = ",")
    xm$Beta[w[1]] <- paste0(unique(xm$Beta[w]), collapse = ",")
    xm$genes[w[1]] <- unique(xm$genes[w]) %>%
      strsplit(., ",") %>%
      unlist() %>%
      unique() %>%
      paste(., collapse = ",")
    xm$Novel[w[1]] <- any(xm$Novel[w])
    xm <- xm[-w[2]]
  }
  xm
}

iga <- read_and_process_ig(snakemake@input$iga)
igg <- read_and_process_ig(snakemake@input$igg)
igm <- read_and_process_ig(snakemake@input$igm)
ig_coloc <- fread(snakemake@input$ig_coloc)

res <- process_ig_and_coloc_results(iga, igg, igm, ig_coloc)

# Remove dubious colocalisation with discordant filtered/unfiltered result
res <- res[rsID != "rs3803800"]

fwrite(res, sep = '\t', file = snakemake@output[[1]])
