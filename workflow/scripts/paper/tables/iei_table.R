library(data.table)

iga <- fread(snakemake@input$iga)
igg <- fread(snakemake@input$igg)
igm <- fread(snakemake@input$igm)
iei <- fread(snakemake@input$iei, select = c('Disease', 'Genetic defect', 'B cell count', 'T cell count', 'Immunoglobulin levels', 'Associated features', 'sanitised_gene', 'OMIM'))

iga[, Isotype := 'IgA']
igm[, Isotype := 'IgM']
igg[, Isotype := 'IgG']

iga[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]
igm[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]
igg[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]

dat <- rbindlist(lapply(list(iga, igm, igg), function(x) x[, .(Isotype, rsID, chromosome, base_pair_location, topGene, iei_gene, distance_to_iei_gene)]))

merged <- merge(dat, iei[, .(Disease, `Genetic defect`, sanitised_gene, `B cell count`, `T cell count`, `Immunoglobulin levels`, `Associated features`, OMIM)], by.x = 'iei_gene', by.y = 'sanitised_gene', all.x = T)

merged <- unique(merged, by = c('Isotype', 'rsID', 'Genetic defect'))

merged <- merged[!(iei_gene %in% c('IGHM', 'IGKC'))]

fwrite(merged[order(Isotype, chromosome, base_pair_location), .(Isotype, rsID, Chromosome = chromosome, Position = base_pair_location, `IEI gene` = `Genetic defect`, `Distance to gene` = distance_to_iei_gene, IEI = Disease, OMIM)], sep = '\t', file = snakemake@output[[1]])

