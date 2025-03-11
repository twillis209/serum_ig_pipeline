library(data.table)
library(MendelianRandomization)
setDTthreads(snakemake@threads)

beta_x <- sprintf("beta.%s", snakemake@isotype)
se_x <- sprintf("se.%s", snakemake@isotype)
beta_y <- sprintf("beta.%s", snakemake@non_ig)
se_y <- sprintf("se.%s", snakemake@non_ig)

merged <- fread(snakemake@input$merged, select = c('rsid', "chromosome", "base_pair_location", beta_x, se_x, beta_y, se_y))

lead_snps <- fread(snakemake@input$lead_snps)

dat <- merge(lead_snps, merged, by = 'rsid', all.x = T)

dat <- dat[!(rsid %in% snakemake@params$snps_to_exclude)]

dat <- na.omit(dat, c(beta_x, beta_y))

mr_input_obj <- mr_input(bx = dat[, beta_x],
                         bxse = dat[, se_x],
                         by = dat[, beta_y],
                         byse = dat[, se_y],
                         snps = dat[, rsid])

res <- mr_allmethods(mr_input_obj, method = 'all')

fwrite(res@Values, file = snakemake@output$tsv, sep = '\t')

png(snakemake@output$png, units = 'in', width = 6, height = 4, res = 300)
mr_plot(res)
dev.off()

png(snakemake@output$png_ivw, units = 'in', width = 6, height = 4, res = 300)
mr_plot(mr_input_obj, line = 'ivw', interactive = F)
dev.off()
