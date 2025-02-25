library(data.table)
library(biomaRt)

dat <- fread(snakemake@input[[1]])

dat[, HPO := NULL]

# Strip whitespace
dat[, sanitised_gene := gsub(' +$', '', sanitised_gene)]
dat[, sanitised_gene := gsub('^ +', '', sanitised_gene)]

# Handle awkward cases
dat[`Genetic defect` %like% "CD40 \\(TNFRSF5\\)", sanitised_gene := 'CD40']
dat[`Genetic defect` == "CD40LG (TNFSF5)", sanitised_gene := 'CD40LG']
dat[`Genetic defect` %like% 'KMT2D', sanitised_gene := 'KMT2D']
dat[`Genetic defect` %like% 'MOGS', sanitised_gene := 'MOGS']
dat[`Genetic defect` %like% 'PIK3CD', sanitised_gene := 'PIK3CD']
dat[`Genetic defect` %like% 'IKBKG', sanitised_gene := 'IKBKG']
dat[`Genetic defect` %like% 'C4A\\+C4B', sanitised_gene := 'C4A']
dat[`Genetic defect` %like% 'CFHR1', sanitised_gene := 'CFHR1']
dat[`Genetic defect` %like% 'IL12B ', sanitised_gene := 'IL12B']

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

query_genes <- function(x) getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = x,
  mart = ensembl
  )

query_genes_by_ensembl_id <- function(x) getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
    filters = "ensembl_gene_id",
    values = x,
    mart = ensembl
    )

res_dat <- data.table(query_genes(dat[!is.na(sanitised_gene), sanitised_gene]))

merged <- merge(dat[!is.na(sanitised_gene)], res_dat, by.x = 'sanitised_gene', by.y = 'hgnc_symbol', all.x = T)

collapsed <- merged[, .(disorder = paste(unique(Disease), collapse = ', ')), by = .(sanitised_gene, ensembl_gene_id, chromosome_name, start_position, end_position)]


# non_hgnc_symbols
list("ADAR1" = "ENSG00000160710",
     "CD20" = "ENSG00000156738",
     "CD21" = "ENSG00000117322",
     "G6PT1" = "ENSG00000137700",
     "IL12B" = "ENSG00000113302",
     "NCKAPIL" =  "ENSG00000123338",
     "NOLA2" = "ENSG00000145912",
     "NOLA3" = "ENSG00000182117",
     "POLE1" = "ENSG00000177084",
     "PSEN" = "ENSG00000080815",
     "SKIV2L" = "ENSG00000204351",
     "TAZ" = "ENSG00000018408",
      "TNFRSF6" = "ENSG00000026103",
      "TNFSF6" = "ENSG00000117560",
     "TTC37" = "ENSG00000198677",
      "XRCC9" = "ENSG00000221829"),


collapsed[, chromosome_name := gsub('HSCHR', '', chromosome_name)]
collapsed[, chromosome_name := gsub('_\\d_CTG\\d+', '', chromosome_name)]
collapsed[, chromosome_name := gsub('6_MHC.+', '6', chromosome_name)]
collapsed[chromosome_name == 'HG1398_PATCH', chromosome_name := '12']
collapsed[chromosome_name == 'HG2217_PATCH', chromosome_name := '11']
collapsed[chromosome_name == 'HG109_PATCH', chromosome_name := '19']
collapsed[chromosome_name == "HG1_PATCH", chromosome_name := '14']
collapsed[chromosome_name == "HG2047_PATCH", chromosome_name := '12']
collapsed[chromosome_name == "HG2526_HG2573_PATCH", chromosome_name := '14']
collapsed[chromosome_name == "HG2334_PATCH", chromosome_name := '10']
collapsed[chromosome_name == "HG1395_PATCH", chromosome_name := '5']
collapsed[chromosome_name == "HG2275_PATCH", chromosome_name := '2']
collapsed[chromosome_name == "21_1", chromosome_name := '21']

fwrite(dat, file = snakemake@output[[1]], sep = '\t')
