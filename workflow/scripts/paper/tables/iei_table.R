library(data.table)

iga <- fread(snakemake@input$iga)
igg <- fread(snakemake@input$igg)
igm <- fread(snakemake@input$igm)

iga[, Isotype := 'IgA']
igm[, Isotype := 'IgM']
igg[, Isotype := 'IgG']

iga[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]
igm[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]
igg[, rsID := paste0(rsid, ':', other_allele, '>', effect_allele)]

dat <- rbindlist(lapply(list(iga, igm, igg), function(x) x[, .(Isotype, rsID, chromosome, base_pair_location, Novel, genes, iei_hgnc_symbol, distance_to_iei_gene, IEIs)]))

dat <- dat[!(iei_hgnc_symbol %in% c('IGHM', 'IGKC'))]

fwrite(dat[order(Isotype, chromosome, base_pair_location), .(Isotype, rsID, Chromosome = chromosome, Position = base_pair_location, Novel, `IEI gene` = iei_hgnc_symbol, `Distance to gene` = distance_to_iei_gene, IEIs)], sep = '\t', file = snakemake@output[[1]])

