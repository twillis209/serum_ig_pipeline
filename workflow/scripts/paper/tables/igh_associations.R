library(data.table)
library(stringr)

consequences <- data.table(mostSevereConsequence = c('missense_variant', 'upstream_gene_variant', 'intron_variant', 'downstream_gene_variant', 'intergenic_variant'), variant_effect = c('missense', 'upstream', 'intronic', 'downstream', 'intergenic'))

dat <- fread(snakemake@input[[1]], sep = '\t', header = T, select = c('dataset', 'rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'beta', 'standard_error', 'p_value', 'nearestGene', 'mostSevereConsequence'))

dat[dataset %like% '-', c('study', 'isotype') := tstrsplit(dataset, split = '-')]
dat[dataset %like% '-', study := snakemake@config$gwas_datasets[[dataset]][['pretty_study']], by = 1:nrow(dat)]
dat[!(dataset %like% '-'), `:=` (study = 'Meta-analysis', isotype = dataset)]
dat[, isotype := snakemake@config$pretty_isotypes[isotype], by = 1:nrow(dat)]

dat[mostSevereConsequence == '', mostSevereConsequence := 'intergenic_variant']

dat <- merge(dat, consequences, all.x = T, by = 'mostSevereConsequence')

dat[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]

dat[, `:=` (Study = study,
            Isotype = isotype,
            rsID = rsID,
            Chromosome = chromosome,
            Position = base_pair_location,
            Beta = beta,
            `Standard error` = standard_error,
            `p-value` = p_value,
            `Nearest gene` = nearestGene,
            `Variant effect` = variant_effect)]

fwrite(dat[order(Isotype, Study), .(Study, Isotype, rsID, Chromosome, Position, Beta, `Standard error`, `p-value`, `Nearest gene`, `Variant effect`)], sep = '\t', file = snakemake@output[[1]])
