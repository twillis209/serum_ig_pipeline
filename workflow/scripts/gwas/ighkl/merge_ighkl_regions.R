library(data.table)
setDTthreads(snakemake@threads)
library(stringr)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
rsid_col <- snakemake@config$rsid_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col

coord_cols <- c(rsid_col, chr_col, bp_col, ref_col, alt_col)
stat_cols <- c(beta_col, se_col, p_col)

studies <- names(snakemake@input)[names(snakemake@input) != ""]

dats <- list()

for (x in studies) {
  dat <- fread(snakemake@input[[x]])
  dat <- dat[((chromosome == snakemake@config$loci[["igh"]]$chrom &
    base_pair_location %between% c(
      snakemake@config$loci[["igh"]]$start - snakemake@params$flank,
      snakemake@config$loci[["igh"]]$stop + snakemake@params$flank
    ))) |
    ((chromosome == snakemake@config$loci[["igk"]]$chrom &
      base_pair_location %between% c(
        snakemake@config$loci[["igk"]]$start - snakemake@params$flank,
        snakemake@config$loci[["igk"]]$stop + snakemake@params$flank
      ))) |
    ((chromosome == snakemake@config$loci[["igl"]]$chrom &
      base_pair_location %between% c(
        snakemake@config$loci[["igl"]]$start - snakemake@params$flank,
        snakemake@config$loci[["igl"]]$stop + snakemake@params$flank
      )))]

      cols_to_append_isotype <- str_subset(names(dat), regex(paste(stat_cols, collapse = "|"), ignore_case = TRUE))

      new_stat_cols <- paste(cols_to_append_isotype, x, sep = ".")

      setnames(dat, cols_to_append_isotype, new_stat_cols)

      cols <- c(coord_cols, new_stat_cols)

      dats[[x]] <- dat[, ..cols]
}

merged <- Reduce(function(x, y) merge(x, y, by = coord_cols, all = TRUE), dats)

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
