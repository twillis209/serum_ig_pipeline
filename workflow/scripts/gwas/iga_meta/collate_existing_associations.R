library(data.table)

ebi <- fread(snakemake@input$ebi, sep = '\t', select = c('CHR_ID', 'CHR_POS', 'SNPS', 'MAPPED_GENE', 'P-VALUE'))
ebi <- ebi[, .(CHR_ID = snakemake@config$chr_col, CHR_POS = snakemake@config$bp_col, SNPS = snakemake@config$rsid_col, MAPPED_GENE = Genes, `P-VALUE` = snakemake@config$p_col)]

# TODO filter ebi for significance, they're not always GWS
# TODO handle liu as special case
# TODO handle willis as special case

epic <- fread(snakemake@input$epic)
epic[, dataset := 'epic']
pietzner <- fread(snakemake@input$pietzner)
pietzner[, dataset := 'pietzner']
gudjonsson <- fread(snakemake@input$gudjonsson)
gudjonsson[, dataset := 'gudjonsson']
eldjarn <- fread(snakemake@input$eldjarn)
eldjarn[, dataset := 'eldjarn']
scepanovic <- fread(snakemake@input$scepanovic)
scepanovic[, dataset := 'scepanovic']
dennis <- fread(snakemake@input$dennis)
dennis[, dataset := 'dennis']

merged <- rbindlist(list(ebi, liu, willis, epic, pietzner, gudjonsson, eldjarn, scepanovic, dennis), fill = T)
