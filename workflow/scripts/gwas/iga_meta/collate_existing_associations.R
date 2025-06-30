library(data.table)

studies <- names(snakemake@input)[names(snakemake@input) != '']

dats <- list()

ebi <- fread(snakemake@input$ebi, sep = '\t', select = c('CHR_ID', 'CHR_POS', 'SNPS', 'MAPPED_GENE', 'P-VALUE'))
ebi <- ebi[, .(chromosome = CHR_ID, base_pair_location = CHR_POS, rsid = SNPS, Genes = MAPPED_GENE, p_value = `P-VALUE`)]
ebi[, dataset := 'ebi']
ebi <- ebi[p_value <= 5e-8]

dats[['ebi']] <- ebi

if('liu' %in% studies) {
  liu <- fread(snakemake@input$liu)
  setnames(liu, 'Locus', 'Genes')
  liu[, dataset := 'liu']
  dats[['liu']] <- liu
}

if('willis' %in% studies) {
  willis <- fread(snakemake@input$willis, select = c('CHR38', 'BP38', 'REF', 'ALT', 'rsID', 'topGene', 'P.meta'))
  willis <- willis[, .(chromosome = CHR38, base_pair_location = BP38, effect_allele = ALT, other_allele = REF, rsid = rsID, Genes = topGene, p_value = P.meta)]
  willis[, dataset := 'willis']
  dats[['willis']] <- willis
}

for(x in studies[!(studies %in% c('ebi', 'liu', 'willis'))]) {
  dat <- fread(snakemake@input[[x]], select = c('chromosome', 'base_pair_location', 'rsid', 'effect_allele', 'other_allele', 'p_value', 'genes'), header = T)
  setnames(dat, 'genes', 'Genes')
  dat[, dataset := x]
  dats[[x]] <- dat
}

merged <- rbindlist(dats, fill = T)

fwrite(merged[order(chromosome, base_pair_location)], file = snakemake@output[[1]], sep = '\t')
