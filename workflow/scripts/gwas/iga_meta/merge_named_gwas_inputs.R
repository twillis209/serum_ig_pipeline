library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
pos_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col
rsid_col <- snakemake@config$rsid_col

cols <- c(chr_col, pos_col, ref_col, alt_col, beta_col, se_col, p_col, rsid_col)
coord_cols <- cols[1:4]
cols_to_relabel <- cols[5:8]

dats <- list()

for(i in seq_along(snakemake@input[names(snakemake@input) != ''])) {
  abbrv <- names(snakemake@input)[names(snakemake@input) != ''][i]
  dat <- fread(snakemake@input[[i]], select = cols)

  setnames(dat, cols_to_relabel, paste(cols_to_relabel, abbrv, sep = '.'))

  dats[[i]] <- dat
}

merged <- Reduce(function(x,y) merge(x, y, by = coord_cols, all = T), dats)

rm(dats)

#merged[, rsid := {
#  rsIDs <- unique(na.omit(unlist(.SD)))
#  if(length(rsIDs) == 1) {
#    rsIDs
#  } else if(length(rsIDs > 1)) {
#    paste(rsIDs, collapse = ',')
#  } else {
#    NA
#  }
#}, .SDcols = patterns("^rsid\\.")]

fwrite(merged, file = snakemake@output[[1]], sep = '\t')
