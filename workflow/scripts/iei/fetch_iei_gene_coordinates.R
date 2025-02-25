library(data.table)
library(biomart)

query_gene_with_ensembl <- function(gene_symbol) {
  base_url = "https://rest.ensembl.org"

  query_url <- paste0(base_url, "/lookup/symbol/human/", gene_symbol, "?content-type=application/json")

  response <- GET(query_url, accept("application/json"))

  if (status_code(response) == 200) {
    # Parse the JSON response
    gene_data <- content(response, as = "parsed", type = "application/json")
    # Print gene details
    print(gene_data)
  } else {
    print(paste("Error:", status_code(response)))
  }

  res <- data.table(t(res))
}

library(biomaRt)

# Use Ensembl as the database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Example list of gene symbols
gene_list <- c("BRCA1", "TP53", "EGFR", "MYC", "KRAS")

# Retrieve Ensembl IDs and other gene information
gene_data <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "hgnc_symbol",
  values = 'BRCA1',
  mart = ensembl
)

# Print result
print(gene_data)




