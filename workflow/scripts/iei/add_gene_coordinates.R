library(data.table)
library(biomaRt)

dat <- fread(snakemake@input[[1]])

dat[, HPO := NULL]

# Strip whitespace
dat[, sanitised_gene := gsub(' +$', '', `Genetic defect`)]
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
dat[`Genetic defect` %like% 'IL12B', sanitised_gene := 'IL12B']

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

# Check res_dat[duplicated(sanitised_gene)], we get multiple ENSGs and coordinates for some symbols, many 
res_dat <- data.table(query_genes(dat[!is.na(sanitised_gene), sanitised_gene]))
res_dat[, sanitised_gene := hgnc_symbol]
uniq_res_dat <- unique(res_dat, by = 'sanitised_gene')

merged <- merge(dat, uniq_res_dat, by = 'sanitised_gene', all.x = T)

# non_hgnc_symbols
genes_without_hgnc <- list(
  "ADAR1" = "ENSG00000160710",
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
  "XRCC9" = "ENSG00000221829")

without_hgnc_dat <- data.table(sanitised_gene = names(genes_without_hgnc), ensembl_gene_id = unlist(genes_without_hgnc))

without_hgnc_res_dat <- data.table(query_genes_by_ensembl_id(without_hgnc_dat[, ensembl_gene_id]))

without_hgnc_res_dat <- without_hgnc_res_dat[without_hgnc_dat, on = 'ensembl_gene_id']

# sanitised gene
merged[without_hgnc_res_dat, on = 'sanitised_gene', `:=` (
                                                      hgnc_symbol = i.hgnc_symbol,
                                                      chromosome_name = i.chromosome_name,
                                                       start_position = i.start_position,
                                                       end_position = i.end_position,
                                                       ensembl_gene_id = i.ensembl_gene_id
                                                    )]

merged[, chromosome_name := gsub('HSCHR', '', chromosome_name)]
merged[, chromosome_name := gsub('_\\d_CTG\\d+', '', chromosome_name)]
merged[, chromosome_name := gsub('6_MHC.+', '6', chromosome_name)]
merged[chromosome_name == 'HG1398_PATCH', chromosome_name := '12']
merged[chromosome_name == 'HG2217_PATCH', chromosome_name := '11']
merged[chromosome_name == 'HG109_PATCH', chromosome_name := '19']
merged[chromosome_name == "HG1_PATCH", chromosome_name := '14']
merged[chromosome_name == "HG2047_PATCH", chromosome_name := '12']
merged[chromosome_name == "HG2526_HG2573_PATCH", chromosome_name := '14']
merged[chromosome_name == "HG2334_PATCH", chromosome_name := '10']
merged[chromosome_name == "HG1395_PATCH", chromosome_name := '5']
merged[chromosome_name == "HG2275_PATCH", chromosome_name := '2']
merged[chromosome_name == "21_1", chromosome_name := '21']

collapsed <- merged[, .(disorder = paste(unique(Disease), collapse = ', ')), by = .(sanitised_gene, hgnc_symbol, ensembl_gene_id, chromosome_name, start_position, end_position)]

fwrite(merged, file = snakemake@output[['sanitised']], sep = '\t')

fwrite(collapsed, file = snakemake@output[['collapsed']], sep = '\t')

