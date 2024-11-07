library(data.table)
setDTthreads(snakemake@threads)
library(EnsDb.Hsapiens.v86)
library(locuszoomr)
library(ggplot2)
library(gridExtra)

dat <- fread(snakemake@input[[1]])

dat[, SNPID := paste(CHR38, BP38, REF, ALT, sep = ':')]

setnames(dat, c('CHR38', 'BP38'), c('chrom', 'pos'))

p_cols <- c('P.liu_decode', 'P.dennis', 'P.lyons', 'P.meta')
subtitles <- c('Liu', 'Dennis', 'Lyons', 'Meta-analysis')

xrange <- c(snakemake@params$index_snp_pos - snakemake@params$flank, snakemake@params$index_snp_pos + snakemake@params$flank)

loci <- lapply(p_cols, function(x) locus(data = dat, p = x, labs = 'SNPID', seqname = snakemake@params$index_snp_seqname, xrange = xrange, ens_db = 'EnsDb.Hsapiens.v86'))

names(loci) <- subtitles

ymax <- max(sapply(loci, function(x) max(x$data$logP, na.rm = T)))

gene_track_pl <- gg_genetracks(loci[[1]])

pls <- lapply(seq_along(loci), function(i) gg_scatter(loci[[i]])+ggtitle(subtitles[i])+ylim(0, ymax)+geom_hline(yintercept = -log10(5e-8), linetype = 'dashed'))

out_pls <- list()

layout_mat <- cbind(1:3, 4:6)
out_pls[1:2] <- pls[1:2]
out_pls[4:5] <- pls[3:4]
out_pls[[3]] <- gene_track_pl
out_pls[[6]] <- gene_track_pl

title_chr <- loci[[1]]$seqname
title_range_left <- format(loci[[1]]$xrange[1], big.mark = ',')
title_range_right <- format(loci[[1]]$xrange[2], big.mark = ',')

ggsave(arrangeGrob(grobs = out_pls, layout_matrix = layout_mat, top = sprintf('%s; %s:%s-%s', snakemake@wildcards[['locus']], title_chr, title_range_left, title_range_right)), file = snakemake@output[[1]], width = 10, height = 12, unit = 'in')
