library(data.table)
setDTthreads(snakemake@threads)
library(stringr)

lead <- fread(snakemake@input$lead)

chr_col <- snakemake@config$chr_col
bp_col <- snakemake@config$bp_col
ref_col <- snakemake@config$ref_col
alt_col <- snakemake@config$alt_col
beta_col <- snakemake@config$beta_col
se_col <- snakemake@config$se_col
p_col <- snakemake@config$p_col

coord_cols <- c(chr_col, bp_col, ref_col, alt_col)

inclusion_wildcards <- snakemake@wildcards[str_detect(names(snakemake@wildcards), '_inclusion')]

studies <- character()

for(x in inclusion_wildcards) {
  if(str_detect(x, 'with_')) {
    studies <- c(studies, (str_split_1(x, '_')[2]))
  }
}

cols <- as.vector(outer(c(beta_col, se_col, p_col), studies, FUN = paste, sep = '.'))

cols <- c(coord_cols, cols)

merged <- fread(snakemake@input$merged, select = cols)

lead_with_sumstats <- merge(lead, merged, by = coord_cols, all.x = T)

fwrite(lead_with_sumstats, file = snakemake@output[[1]], sep = '\t')
