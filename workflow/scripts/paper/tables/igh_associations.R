library(data.table)
library(stringr)

consequences <- data.table(most_severe_consequence = c('missense_variant', 'upstream_gene_variant', 'intron_variant', 'downstream_gene_variant', 'intergenic_variant'), variant_effect = c('missense', 'upstream', 'intronic', 'downstream', 'intergenic'))

igh <- fread(snakemake@input[['igh']], sep = '\t', header = T, select = c('dataset', 'rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'beta', 'standard_error', 'p_value', 'nearest_gene', 'most_severe_consequence'))
igk <- fread(snakemake@input[['igk']], sep = '\t', header = T, select = c('dataset', 'rsid', 'chromosome', 'base_pair_location', 'other_allele', 'effect_allele', 'beta', 'standard_error', 'p_value', 'nearest_gene', 'most_severe_consequence'))

dat <- rbindlist(list(igh, igk), fill = T)

dat[dataset %like% '-', c('study', 'isotype') := tstrsplit(dataset, split = '-')]

pretty_study <- sapply(dat[dataset %like% '-', unique(dataset)], function(x) snakemake@config$gwas_datasets[[x]][['pretty_study']])
pretty_study_dat <- data.table(study = names(pretty_study), pretty_study = pretty_study)

dat <- merge(dat, pretty_study_dat, all.x = T, by.x = 'dataset', by.y = 'study')

dat[!(dataset %like% '-'), `:=` (pretty_study = 'Meta-analysis', isotype = dataset)]
dat[, isotype := snakemake@config$pretty_isotypes[isotype], by = 1:nrow(dat)]

dat[most_severe_consequence == '', most_severe_consequence := 'intergenic_variant']

dat <- merge(dat, consequences, all.x = T, by = 'most_severe_consequence')

dat[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]

dat[, `:=` (Study = pretty_study,
            Isotype = isotype,
            Variant = rsID,
            Chromosome = chromosome,
            Position = base_pair_location,
            Beta = beta,
            `Standard error` = standard_error,
            `p-value` = p_value,
            `Nearest gene` = nearest_gene,
            `Variant effect` = variant_effect)]

fwrite(dat[order(Isotype, Study), .(Study, Isotype, rsID, Chromosome, Position, Beta, `Standard error`, `p-value`, `Nearest gene`, `Variant effect`)], sep = '\t', file = snakemake@output[[1]])
