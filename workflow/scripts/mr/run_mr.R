library(data.table)
library(MendelianRandomization)
setDTthreads(snakemake@threads)

beta_y <- sprintf("beta.%s-meta", snakemake@wildcards$isotype)
se_y <- sprintf("standard_error.%s-meta", snakemake@wildcards$isotype)
beta_x <- sprintf("beta.%s", snakemake@wildcards$non_ig)
se_x <- sprintf("standard_error.%s", snakemake@wildcards$non_ig)

merged <- fread(snakemake@input$merged, select = c('rsid', "chromosome", "base_pair_location", beta_x, se_x, beta_y, se_y))

lead_snps <- fread(snakemake@input$instruments)

dat <- merge(lead_snps, merged, by = 'rsid')

dat <- dat[!(rsid %in% snakemake@params$snps_to_exclude)]

dat <- na.omit(dat, c(beta_x, beta_y))

mr_input_obj <- mr_input(bx = dat[, get(beta_x)],
                         bxse = dat[, get(se_x)],
                         by = dat[, get(beta_y)],
                         byse = dat[, get(se_y)],
                         snps = dat[, rsid])

saveRDS(mr_input_obj, snakemake@output$mr_input)

res <- mr_allmethods(mr_input_obj, method = 'all')

fwrite(res@Values, file = snakemake@output$tsv, sep = '\t')

png(snakemake@output$png, units = 'in', width = 6, height = 4, res = 300)
mr_plot(res)
dev.off()

png(snakemake@output$png_ivw, units = 'in', width = 6, height = 4, res = 300)
mr_plot(mr_input_obj, line = 'ivw', interactive = F)
dev.off()
