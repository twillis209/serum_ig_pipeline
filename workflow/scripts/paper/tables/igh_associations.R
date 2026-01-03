library(data.table)
library(stringr)
library(GenomicRanges)
library(ensembldb)

edb = EnsDb(snakemake@input[['edb']])

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

## dat[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]

# Replace old nearest genes with Ensembl 113 ones
snps_gr <- GRanges(
  seqnames = dat$chromosome,
  ranges = IRanges(start = dat$base_pair_location, width = 1),
  rsid = dat$rsid
)

genes_gr <- genes(edb)

suppressWarnings(nearest_hits <- distanceToNearest(snps_gr, genes_gr))

nearest_genes <- data.table(
  Variant = mcols(snps_gr[queryHits(nearest_hits)])$rsid,
  nearest_gene_id        = names(genes_gr[subjectHits(nearest_hits)]),
  nearest_gene_name      = genes_gr[subjectHits(nearest_hits)]$gene_name,
  distance_bp    = mcols(nearest_hits)$distance
)

dat[nearest_genes, on = 'Variant', `:=` (nearest_gene = nearest_gene_name, distance_bp = distance_bp)]

dat[, `:=`(
  Study = pretty_study,
  Isotype = isotype,
  rsID = rsid,
  Chromosome = chromosome,
  Position = base_pair_location,
  `Effect allele` = effect_allele,
  `Other allele` = other_allele,
  Beta = beta,
  `Standard error` = standard_error,
  `p-value` = p_value,
  `Nearest gene` = nearest_gene,
  `Distance to nearest gene` = distance_bp,
  `Variant effect` = variant_effect
)]

# Drop dubious variants still in LD with lead SNP
# IgA IGH
dat <- dat[rsID != "rs115597735"]

fwrite(dat[order(Isotype, Study), .(Study, Isotype, rsID, Chromosome, Position, `Effect allele`, `Other allele`, Beta, `Standard error`, `p-value`, `Nearest gene`, `Variant effect`)], sep = '\t', file = snakemake@output[[1]])
