library(data.table)
setDTthreads(snakemake@threads)
library(serumIgPipelineCode)
library(stringr)

dat <- fread(snakemake@input[[1]])

add_study <- function(study_name) {
  update_meta_analysis_estimates(dat, study_name,
                                 effect_label = snakemake@config$beta_col,
                                 se_label = snakemake@config$se_col,
                                 p_label = snakemake@config$p_col)
}

inclusion_wildcards <- snakemake@wildcards[str_detect(names(snakemake@wildcards), '_inclusion')]

for(x in inclusion_wildcards) {
  if(str_detect(x, 'with_')) {
    add_study(str_split_1(x, '_')[2])
  }
}

cols <- c(snakemake@config$chr_col,
          snakemake@config$bp_col,
          snakemake@config$ref_col,
          snakemake@config$alt_col,
          snakemake@config$beta_col,
          snakemake@config$se_col,
          snakemake@config$p_col
          )

fwrite(dat[, ..cols], file = snakemake@output[[1]], sep = '\t')
