library(data.table)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            axis.text = element_text(size = 10, color = "black"),
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10)
          )
          )

dat <- fread(snakemake@input[[1]])

dat[, trim_suffix := ifelse(trimmed, 'trim', 'raw')]

cast_dat <- dcast(dat, first_trait + second_trait + first_snp + second_snp + genes.first_snp + genes.second_snp + r2 ~ trim_suffix, value.var = c('PP.H0.abf', 'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf'), sep = '.')

cast_outlier_dat <- cast_dat[(first_snp == 'rs3803800' & second_snp == 'rs758641530') |
                             (first_snp == 'rs1948915' & second_snp == 'rs747624575')]

cast_outlier_dat[, label := paste(first_snp, second_snp, sep = '-')]

ggsave(ggplot(cast_dat)+
       geom_point(aes(x = PP.H4.abf.raw, y = PP.H4.abf.trim))+
       geom_label_repel(data = cast_outlier_dat, aes(label = label, x = PP.H4.abf.raw, y = PP.H4.abf.trim), color = 'blue', nudge_x = 0.1),
       file = snakemake@output[['h4']])
