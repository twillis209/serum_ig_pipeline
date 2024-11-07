library(data.table)
#library(pidPipelineCode)

# NB: K is prevalence, P is proportion of cases in the study
obs_scale_to_liab_scale <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}

rg <- fread(snakemake@input[['rg']])

metadat <- fread(snakemake@input[['metadata']], select = c('abbrv', 'prevalence', 'case_prop'))

metadat <- metadat[abbrv == snakemake@wildcards[['trait']]]

rg[, `:=` (h2.A.liab.sr = obs_scale_to_liab_scale(h2.A.obs.sr, P = metadat$case_prop, K = metadat$prevalence),
           h2.A.liab.se = obs_scale_to_liab_scale(h2.A.obs.se.sr, P = metadat$case_prop, K = metadat$prevalence))]

fwrite(rg, file = snakemake@output[[1]], sep = '\t')
