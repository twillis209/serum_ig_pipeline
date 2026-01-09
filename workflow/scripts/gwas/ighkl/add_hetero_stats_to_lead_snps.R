library(data.table)

calc_Q_I2 <- function(row) {
  isotype <- row[, tolower(Isotype)]

  config_key <- paste(isotype, "studies", sep = "_")
  studies <- snakemake@config[[config_key]]
  beta_cols <- paste("beta", studies, isotype, sep = ".")
  se_cols <- paste("standard_error", studies, isotype, sep = ".")

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
    p_value <- pchisq(Q, df, lower.tail = FALSE)
    I2 <- max(0, (Q - df) / Q) * 100

    c(Q = Q, df = df, I2 = I2, Q.p_value = p_value)
  }
}

lead_snps <- fread(snakemake@input[['lead_snps']], sep = '\t')
sumstats <- fread(snakemake@input[['sumstats']], sep = '\t')

merged <- merge(lead_snps, sumstats, by = "rsid", all.x = T)

merged[Study == "Meta-analysis", c("Q", "df", "I2", "Q.p_value") :=
  as.list(
    as.data.table(
      t(vapply(.I, function(i) calc_Q_I2(lead_snps[i]), numeric(4)))
    )
  )]

cols_to_keep <- c(names(lead_snps), "Q", "df", "I2", "Q.p_value")

fwrite(merged[, ..cols_to_keep], file = snakemake@output[[1]], sep = "\t")
