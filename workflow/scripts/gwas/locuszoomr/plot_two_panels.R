library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(gridExtra)
library(stringr)

snps_to_label <- c(snakemake@params$snps_to_label_for_A, snakemake@params$snps_to_label_for_B)
snps_to_label <- snps_to_label[!duplicated(snps_to_label)]

relabel_snps <- function(dat, snps_to_label) {
  setkey(dat, SNPID)

  snp_map <- data.table(old = as.character(names(snps_to_label)),
                        new = as.character(snps_to_label))
  setkey(snp_map, old)

  dat[snp_map, SNPID := new]
}

gwas_A  <- fread(snakemake@input[['gwas_A']])
relabel_snps(gwas_A, snps_to_label)
gwas_B  <- fread(snakemake@input[['gwas_B']])
relabel_snps(gwas_B, snps_to_label)

gwas_A <- na.omit(gwas_A, snakemake@params$p_col_A)
gwas_B <- na.omit(gwas_B, snakemake@params$p_col_B)

loc_A <- locus(data = gwas_A, chrom = 'CHR38', pos = 'BP38', p = snakemake@params$p_col_A, ens_db = "EnsDb.Hsapiens.v86", seqname = as.integer(snakemake@params$chrom), xrange = snakemake@params$xrange, index_snp = snps_to_label, labs = 'SNPID')
loc_B <- locus(data = gwas_B, chrom = 'CHR38', pos = 'BP38', p = snakemake@params$p_col_B, ens_db = "EnsDb.Hsapiens.v86", seqname = as.integer(snakemake@params$chrom), xrange = snakemake@params$xrange, index_snp = snps_to_label, labs = 'SNPID')

pls <- list()

pls[[1]] <- gg_scatter(loc_A, labels = snps_to_label, showLD = F, min.segment.length = 0, nudge_y = 5)+ggtitle(snakemake@params$pretty_A)
pls[[2]] <- gg_scatter(loc_B, labels = snps_to_label, showLD = F, min.segment.length = 0, nudge_y = 5)+ggtitle(snakemake@params$pretty_B)
pls[[3]] <- gg_genetracks(loc_A, filter_gene_biotype = snakemake@params$filter_gene_biotype)

ggsave(arrangeGrob(grobs = pls, layout_matrix = as.matrix(1:3),
                   top = sprintf('%s; %s:%s-%s',
                                 str_to_upper(snakemake@wildcards$locus),
                                 snakemake@params$chrom,
                                 format(snakemake@params$xrange[1], big.mark = ',', scientific = F),
                                 format(snakemake@params$xrange[2], big.mark = ',', scientific = F)
                                 )
                   ), file = snakemake@output[[1]], width = 6, height = 9, unit = 'in')
