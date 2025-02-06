library(data.table)
library(stringr)

studies <- names(snakemake@input)[names(snakemake@input) != '']

pretty_isotypes = list(iga = 'IgA', igm = 'IgM', igg = 'IgG')

consequences = list(missense_variant = 'missense', upstream_gene_variant = 'upstream', intronic_variant = 'intronic')

dats <- list()

for(x in studies) {
  dat <- fread(snakemake@input[[x]], sep = '\t', header = T, select = c('rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'beta', 'standard_error', 'p_value', 'nearestGene', 'mostSevereConsequence'))
  study_split <- str_split_1(x, '_')

  pretty_study <- snakemake@config[['gwas_datasets']][[paste(study_split, collapse = '-')]][['pretty_study']]

  dat[, `Variant effect` := '']

  dat[mostSevereConsequence != '', `Variant effect` := consequences[[mostSevereConsequence]]]
  dat[mostSevereConsequence == '', `Variant effect` := 'intergenic']

  dat[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]

  dat[, `:=` (Study = pretty_study,
              Isotype = pretty_isotypes[study_split[2]],
              rsID = rsID,
              Chromosome = chromosome,
              Position = base_pair_location,
              Beta = beta,
              `Standard error` = standard_error,
              `p-value` = p_value,
              `Nearest gene` = nearestGene,
              `Variant effect` = mostSevereConsequence)]

  dats[[x]] <- dat
}

fwrite(rbindlist(dats)[, .(Study, Isotype, rsID, Chromosome, Position, Beta, `Standard error`, `p-value`, `Nearest gene`, `Variant effect`)], sep = '\t', file = snakemake@output[[1]])
