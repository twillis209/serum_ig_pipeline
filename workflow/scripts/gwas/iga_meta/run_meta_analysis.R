library(data.table)
setDTthreads(snakemake@threads)
library(serumIgPipelineCode)

dat <- fread(snakemake@input[[1]])

add_study <- function(study_name) {
  update_meta_analysis_estimates(dat, study_name,
                                 effect_label = snakemake@config$beta_col,
                                 se_label = snakemake@config$se_col,
                                 p_label = snakemake@config$p_col)
}

if(snakemake@wildcards$epic_inclusion == 'with_epic') {
  add_study('epic')
}
if(snakemake@wildcards$liu_inclusion == 'with_liu') {
  add_study('liu')
}
if(snakemake@wildcards$scepanovic_inclusion == 'with_scepanovic') {
  add_study('scepanovic')
}
if(snakemake@wildcards$dennis_inclusion == 'with_dennis') {
  add_study('dennis')
}
if(snakemake@wildcards$pietzner_inclusion == 'with_pietzner') {
  add_study('pietzner')
}
if(snakemake@wildcards$gudjonsson_inclusion == 'with_gudjonsson') {
  add_study('gudjonsson')
}
if(snakemake@wildcards$eldjarn_inclusion == 'with_eldjarn') {
  add_study('eldjarn')
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
