library(data.table)

lead <- fread(snakemake@input[[1]])

lead[, rsID := paste(paste(rsid, other_allele, sep = ':'), effect_allele, sep = '>')]

lead[nearestGene == '' & topGene != '', gene := topGene]
lead[nearestGene != '' & topGene == '', gene := nearestGene]
lead[nearestGene == topGene, gene := nearestGene]
lead[nearestGene != '' & topGene != '' & nearestGene != topGene, gene := paste(nearestGene, topGene, sep = ', ')]

fwrite(lead[, .(rsID = rsID,
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
