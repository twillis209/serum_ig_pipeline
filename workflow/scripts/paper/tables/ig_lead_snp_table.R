library(data.table)

lead <- fread(snakemake@input[[1]])

## lead[, rsID := paste(paste(rsid, other_allele, sep = ':'), effect_allele, sep = '>')]

lead[, chromosome := as.character(chromosome)]

lead[chromosome == "23", chromosome := "X"]

cohorts <- gsub("beta.", "", names(lead)[names(lead) %like% "beta\\.\\w+"])

beta_cols <- paste0("beta.", cohorts)

se_cols <- paste0("standard_error.", cohorts)

calc_Q_I2 <- function(row) {
  beta_cols <- paste0("beta.", cohorts)
  se_cols <- paste0("standard_error.", cohorts)

  beta <- as.numeric(unlist(row[, ..beta_cols]))
  se <- as.numeric(unlist(row[, ..se_cols]))

  # Remove missing studies
  keep <- !is.na(beta) & !is.na(se)
  beta <- beta[keep]
  se <- se[keep]

  k <- length(beta)

  if (k < 2) {
    c(Q = NA_real_, df = k, I2 = NA_real_)
  } else {
    w <- 1 / se^2
    Q <- sum(w * (beta - row[, beta])^2)

    df <- k - 1
    I2 <- max(0, (Q - df) / Q) * 100

    c(Q = Q, df = df, I2 = I2)
  }
}

lead[, c("Q", "df", "I2") :=
  as.list(
    as.data.table(
      t(vapply(.I, function(i) calc_Q_I2(lead[i]), numeric(3)))
    )
  )]


fwrite(lead[, .(rsID = rsid,
               Chromosome = chromosome,
               Position = base_pair_location,
               `Effect allele` = effect_allele,
                `Other allele` = other_allele,
               MAF = maf.meta,
               `Sample size` = sample_size,
               `Gene(s)` = genes,
               `Nearest gene` = nearest_gene_name,
               `Distance to nearest gene` = distance_bp,
               `Missense gene` = missense_gene,
               `QTL genes` = qtl_genes,
               Novel = Novel,
               Beta = beta,
               `Standard error` = standard_error,
               `p-value` = p_value,
               Q = Q,
               `Degrees of freedom` = df,
               I2 = I2
               )],
       file = snakemake@output[[1]],
       sep = '\t')
