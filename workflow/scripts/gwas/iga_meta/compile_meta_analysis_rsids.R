library(data.table)
setDTthreads(snakemake@threads)
library(stringr)

save.image('rsids.RData')

studies <- names(snakemake@input)[names(snakemake@input) != '']

coord_cols <- c(snakemake@config$chr_col,
                snakemake@config$bp_col,
                snakemake@config$ref_col,
                snakemake@config$alt_col)

cols <- c(coord_cols, snakemake@config$rsid_col)

dats <- list()

inclusion_wildcards <- snakemake@wildcards[str_detect(names(snakemake@wildcards), '_inclusion')]

for(x in inclusion_wildcards) {
  if(str_detect(x, 'with_')) {
    study_name <- str_split_1(x, '_')[2]

    dat <- fread(snakemake@input[[study_name]], sep = '\t', select = cols)

    setnames(dat, snakemake@config$rsid_col, paste(snakemake@config$rsid_col, study_name, sep = '.'))

    dats[[study_name]] <- dat
  }
}

merged <- Reduce(function(x, y) merge(x, y, by = coord_cols, all = TRUE), dats)

rm(dats)
gc()

merged[, rsid := fcoalesce(.SD), .SDcols = patterns('^rsid')]

out_cols <- c(coord_cols, 'rsid')

fwrite(merged[, ..out_cols], sep = '\t', file = snakemake@output[[1]])
