library(data.table)
setDTthreads(snakemake@threads)

library(ggplot2)
theme_set(theme_bw()+
          theme(
            axis.title = element_text(size=12),
            plot.title=element_text(hjust=0.5, size=12),
            strip.text=element_text(size=10),
            axis.text.x=element_text(size=6, angle=90, color="black"),
            axis.text.y=element_text(size=10, color="black"),
            legend.title=element_text(size=10),
            legend.text=element_text(size=10)
          )
          )

stratified_qqplot <- function(dat, prin_value_label, cond_value_label = NULL, thresholds = c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {

  dat[, negLogP := -log10(get(prin_value_label))]

  if(is.null(cond_value_label)) {
    dat <- dat[order(get(prin_value_label))]
    dat[, pp := -log10(ppoints(nrow(dat)))]
    dat[, threshold := factor(c(1))]
  } else {

    dats <- list()

    for(i in seq_along(thresholds)) {
      temp_dat <- dat[get(cond_value_label) < thresholds[i]]
      temp_dat <- temp_dat[order(get(prin_value_label))]
      temp_dat[, pp := -log10(ppoints(nrow(temp_dat)))]
      temp_dat[, threshold := factor(thresholds)[i]]
      dats[[i]] <- temp_dat
    }

    dat <- rbindlist(dats)
  }

  ggplot(data = dat)+
    geom_point(aes(x = pp, y = negLogP, group = threshold, colour = threshold))+
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")
}

if(is.null(snakemake@params[['aux_col']])) {
  dat <- fread(snakemake@input[['gwas']], select = snakemake@params[['prin_col']])
  gif_dat <- fread(snakemake@input[['gif']], sep = '\t', header = T)
  gif <- gif_dat[1, lambda_0_50]
  pl <- stratified_qqplot(dat, prin_value_label = snakemake@params[['prin_col']])+
    theme(legend.position = 'none')+
    ylab('Observed -log10(p)')+
    xlab('Expected -log10(p)')+
    ggtitle(sprintf('GWAS of %s (%s et al.), GIF = %.3f', snakemake@params[['pretty_trait']], snakemake@params[['pretty_author']], gif))+
    ylim(0, 30)+
    xlim(0, 10)
} else {
  dat <- fread(snakemake@input[['gwas']], select = c(snakemake@params[['prin_col']], snakemake@params[['aux_col']]))
  pl <- stratified_qqplot(dat, prin_value_label = snakemake@params[['prin_col']], cond_value_label = snakemake@params[['aux_col']])
}

ggsave(pl, file = snakemake@output[[1]])
