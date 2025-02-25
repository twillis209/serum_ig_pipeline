library(data.table)

dat <- fread(snakemake@input[[1]], strip.white = T)

# Omit unknown and chromosomal aberrations
dat[!(`Genetic defect` %like% 'Unknown|\\d+(p|q)'), gene := `Genetic defect`]

# Strip whitespace
dat[, gene := gsub(' +$', '', gene)]

# Handle awkward cases
dat[`Genetic defect` %like% "CD40 \\(TNFRSF5\\)", gene := 'CD40']
dat[`Genetic defect` == "CD40LG (TNFSF5)", gene := 'CD40LG']
dat[`Genetic defect` %like% 'KMT2D', gene := 'KMT2D']
dat[`Genetic defect` %like% 'MOGS', gene := 'MOGS']
dat[`Genetic defect` %like% 'PIK3CD', gene := 'PIK3CD']
dat[`Genetic defect` %like% 'IKBKG', gene := 'IKBKG']

dat[`Genetic defect` %like% 'C4A\\+C4B', gene := 'C4A']
dat[`Genetic defect` %like% 'CFHR1', gene := 'CFHR1']

setnames(dat, 'gene', 'sanitised_gene')

fwrite(dat[sanitised_gene != '', .(tangye_gene = `Genetic defect`, sanitised_gene)], sep = '\t', file = snakemake@output[['gene_list']])

fwrite(dat, sep = '\t', file = snakemake@output[['sanitised']])
