library(data.table)

lead <- fread(snakemake@input[['lead']])

novel <- fread(snakemake@input[['novel']])

lead[, rsID := paste(paste(rsid, other_allele, sep = ':'), effect_allele, sep = '>')]

novel[, Novel := TRUE]

merged <- merge(lead, novel[, .(chromosome, base_pair_location, Novel)], by = c('chromosome', 'base_pair_location'), all.x = T)

merged[is.na(Novel), Novel := FALSE]

merged[nearestGene == '' & topGene != '', gene := topGene]
merged[nearestGene != '' & topGene == '', gene := nearestGene]
merged[nearestGene == topGene, gene := nearestGene]
merged[nearestGene != '' & topGene != '' & nearestGene != topGene, gene := paste(nearestGene, topGene, sep = ', ')]

# TODO not actually sure which is the minor allele
fwrite(merged[, .(rsID = rsID,
               Chromosome = chromosome,
               Position = base_pair_location,
               MAF = maf.meta,
               `Gene(s)` = gene,
               Novel = Novel,
               Beta = beta,
               `Standard error` = standard_error,
               `p-value` = p_value
               )],
       file = snakemake@output[[1]],
       sep = '\t')
