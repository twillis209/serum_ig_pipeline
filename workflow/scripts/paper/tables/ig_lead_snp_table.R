library(data.table)

lead <- fread(snakemake@input[[1]])

lead[, rsID := paste(paste(rsid, other_allele, sep = ':'), effect_allele, sep = '>')]

fwrite(lead[, .(rsID = rsID,
               Chromosome = chromosome,
               Position = base_pair_location,
               MAF = maf.meta,
               `Gene(s)` = genes,
               `Nearest gene` = nearest_gene_name,
               `Distance to nearest gene` = distance_bp,
               `Missense gene` = missense_gene,
               `QTL genes` = qtl_genes,
               Novel = Novel,
               Beta = beta,
               `Standard error` = standard_error,
               `p-value` = p_value
               )],
       file = snakemake@output[[1]],
       sep = '\t')
