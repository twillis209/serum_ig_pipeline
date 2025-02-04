library(data.table)

ebi <- fread(snakemake@input$ebi, sep = '\t', select = c('CHR_ID', 'CHR_POS', 'SNPS', 'MAPPED_GENE', 'P-VALUE'))
ebi <- ebi[, .(chromosome = CHR_ID, base_pair_location = CHR_POS, rsid = SNPS, Genes = MAPPED_GENE, p_value = `P-VALUE`)]
ebi[, dataset := 'ebi']
ebi <- ebi[p_value <= 5e-8]

dats <- list()

for(x in c('pietzner', 'gudjonsson', 'eldjarn', 'scepanovic', 'dennis')) {
  print(x)
  dat <- fread(snakemake@input[[x]], sep = '\t', header = T, select = c('chromosome', 'base_pair_location', 'rsid', 'effect_allele', 'other_allele', 'p_value', 'topGene'))
  setnames(dat, 'topGene', 'Genes')
  dat[, dataset := x]
  dats[[x]] <- dat
}

merged <- rbindlist(c(list(ebi), dats), fill = T)

fwrite(merged[order(chromosome, base_pair_location)], file = snakemake@output[[1]], sep = '\t')
