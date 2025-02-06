library(data.table)
setDTthreads(snakemake@threads)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
rsid_col <- snakemake@config$rsid_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col

cols <- c(chr_col, bp_col, ref_col, alt_col, rsid_col, beta_col, se_col, p_col)
coord_cols <- cols[1:4]
cols_to_relabel <- cols[5:8]

dats <- list()

studies <- names(snakemake@input)[names(snakemake@input) != ""]

for(x in studies) {
  dat <- fread(snakemake@input[[x]], select = cols)

  setnames(dat, cols_to_relabel, paste(cols_to_relabel, x, sep = "."))

  dats[[x]] <- dat
}

if (is.null(snakemake@wildcards$join)) {
  merged <- Reduce(function(x, y) merge(x, y, by = coord_cols, all = TRUE), dats)
} else if (snakemake@wildcards$join == "outer") {
  merged <- Reduce(function(x, y) merge(x, y, by = coord_cols, all = TRUE), dats)
} else if (snakemake@wildcards$join == "inner") {
  merged <- Reduce(function(x, y) merge(x, y, by = coord_cols), dats)
} else {
  stop("Invalid join type")
}

if (!is.null(snakemake@wildcards$variant_set)) {
  if (snakemake@wildcards$variant_set == "sans_mhc") {
    merged <- merged[!(as.integer(chr) == 6 & bp %between% c(24e6, 45e6)),
                     env = list(chr = chr_col, bp = bp_col)]
  }
}

merged[, rsid := fcoalesce(.SD), .SDcols = patterns('^rsid')]

out_cols <- !(names(merged) %like% '^rsid\\.')

fwrite(merged[, ..out_cols], file = snakemake@output[[1]], sep = "\t")
