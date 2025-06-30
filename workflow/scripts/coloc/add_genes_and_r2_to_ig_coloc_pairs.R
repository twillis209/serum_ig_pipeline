library(data.table)
library(LDlinkR)

all_pairs <- fread(snakemake@input$all_pairs)

iga <- fread(snakemake@input$iga)
igg <- fread(snakemake@input$igg)
igm <- fread(snakemake@input$igm)

rbound_genes <- rbindlist(list(iga[, .(rsid, genes)],
                               igg[, .(rsid, genes)],
                               igm[, .(rsid, genes)]))

rbound_genes <- rbound_genes[, .(genes = paste(unique(unlist(strsplit(genes, ","))), collapse = ",")), by = rsid][genes != '']

setkey(rbound_genes, rsid)

setkey(all_pairs, first_snp)
all_pairs[rbound_genes, genes.first_snp := genes]

setkey(all_pairs, second_snp)
all_pairs[rbound_genes, genes.second_snp := genes]

all_pairs[, max_post := names(.SD)[max.col(.SD, ties.method = 'first')], .SDcols = patterns('PP.H')]
all_pairs[, max_post := gsub('PP\\.|\\.abf', '', max_post)]

for(i in seq_len(all_pairs[, .N])) {
  r2 <- tryCatch({
  ld <- LDpop(
    var1 = all_pairs[i, first_snp],
    var2 = all_pairs[i, second_snp],
    pop = "EUR",
    r2d = "r2",
    token = snakemake@config$LDlink$token,
    file = FALSE,
    genome_build = "grch38_high_coverage"
  )
  ld[ld$Abbrev == 'EUR', 'R2']
  },
  error = function(e) { message(e); NA })

 set(all_pairs, i, "r2", r2)
}

fwrite(all_pairs, file = snakemake@output[[1]], sep = '\t')
