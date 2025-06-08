library(gwasrapidd)
library(dplyr)
library(tidyr)
library(fuzzyjoin)

# Time-consuming function btw, especially when looking for lymphocyte count associations
get_gws_associations_for_lymphocyte_trait <- function(reported_trait) {
  reported_traits_studies <- get_studies(reported_trait = reported_trait)

  reported_traits_assocs <- get_associations(study_id = reported_traits_studies@studies$study_id)

  gws_assocs <- filter(reported_traits_assocs@associations, pvalue < 5e-8) %>%
    drop_na(pvalue)

  gws_vars <- get_variants(association_id = gws_assocs$association_id)

  gws_assoc_gws_vars <- gws_assocs %>%
    mutate(
      beta = beta_number * if_else(beta_direction == "decrease", -1, 1)
    ) %>%
    select(
      association_id,
      pvalue,
      beta,
      standard_error
    ) %>%
    left_join(
      reported_traits_assocs@risk_alleles[, c("association_id", "variant_id")],
      by = "association_id"
    ) %>%
    left_join(
      gws_vars@variants[, c("variant_id", "chromosome_name", "chromosome_position")],
      by = "variant_id"
    ) %>%
    left_join(
      group_by(gws_vars@ensembl_ids, variant_id) %>% summarize(genes = paste(gene_name, collapse = ","))
    ) %>%
    rename(
      chromosome = chromosome_name,
      base_pair_location = chromosome_position
    )

  gws_assoc_gws_vars_studies <- association_to_study(gws_assoc_gws_vars$association_id)

  gws_assoc_gws_vars <- left_join(
    gws_assoc_gws_vars,
    gws_assoc_gws_vars_studies,
    by = "association_id"
  ) %>%
    left_join(
      reported_traits_studies@studies[, c("study_id", "reported_trait")],
      by = "study_id"
    )
}

fuzzy_join_assocs <- function(a_assocs, b_assocs, a_suffix, b_suffix, max_distance = 1e6) {
  assocs <- fuzzy_inner_join(
    select(a_assocs, reported_trait, variant_id, chromosome, base_pair_location, genes, beta),
    select(b_assocs, reported_trait, variant_id, chromosome, base_pair_location, genes, beta),
    by = c("chromosome" = "chromosome", "base_pair_location" = "base_pair_location"),
    match_fun = list(`==`, function(x, y) abs(x - y) <= max_distance)
  )

  names(assocs) <- gsub(".x", sprintf(".%s", a_suffix), names(assocs))
  names(assocs) <- gsub(".y", sprintf(".%s", b_suffix), names(assocs))

  assocs
}

leukocyte_count_efo_id <- "OBA_VT0000217"
lymphocyte_count_efo_id <- "EFO_0004587"

leukocyte_count_studies <- get_studies(efo_id = leukocyte_count_efo_id)
lymphocyte_count_studies <- get_studies(efo_id = lymphocyte_count_efo_id)

# CD19-positive B-lymphocyte count is a child trait for leukocyte count but does not appear as a reported trait, also not a child trait for lymphocyte count
leukocyte_count_b_cell_reported_traits <- leukocyte_count_studies@studies %>%
  dplyr::filter(grepl("^B cell", reported_trait, ignore.case = TRUE)) %>%
  dplyr::pull(reported_trait) %>%
  unique()
lymphocyte_count_b_cell_reported_traits <- lymphocyte_count_studies@studies %>%
  dplyr::filter(grepl("B cell", reported_trait, ignore.case = TRUE)) %>%
  dplyr::pull(reported_trait) %>%
  unique()

misc_b_cell_traits <- c("CD19-positive B-lymphocyte count")

b_cell_reported_traits <- c(misc_b_cell_traits, leukocyte_count_b_cell_reported_traits, lymphocyte_count_b_cell_reported_traits)

gws_b_cell_assocs <- get_gws_associations_for_lymphocyte_trait(b_cell_reported_traits)

leukocyte_count_t_cell_reported_traits <- leukocyte_count_studies@studies %>%
  dplyr::filter(grepl("T cell", reported_trait, ignore.case = TRUE)) %>%
  dplyr::pull(reported_trait) %>%
  unique()

lymphocyte_count_t_cell_reported_traits <- lymphocyte_count_studies@studies %>%
  dplyr::filter(grepl("T cell", reported_trait, ignore.case = TRUE)) %>%
  dplyr::pull(reported_trait) %>%
  unique()

misc_t_cell_traits <- c("CD3- lymphocyte Absolute Count")

t_cell_reported_traits <- c(misc_t_cell_traits, leukocyte_count_t_cell_reported_traits, lymphocyte_count_t_cell_reported_traits)

gws_t_cell_assocs <- get_gws_associations_for_lymphocyte_trait(t_cell_reported_traits)

total_lymphocyte_count_traits <- lymphocyte_count_studies@studies %>%
  dplyr::filter(!grepl("T|B cell", reported_trait)) %>%
  dplyr::filter(!grepl("plasma cell", reported_trait, ignore.case = T)) %>%
  dplyr::filter(!grepl("gene-based burden", reported_trait, ignore.case = T)) %>%
  dplyr::filter(!grepl("fraction", reported_trait, ignore.case = T)) %>%
  dplyr::filter(!grepl("variance|CNV", reported_trait, ignore.case = T)) %>%
  dplyr::filter(!grepl("short tandem repeats", reported_trait, ignore.case = T)) %>%
  dplyr::filter(grepl("lymphocyte", reported_trait, ignore.case = T)) %>%
  dplyr::filter(!grepl("percent", reported_trait, ignore.case = T)) %>%
  pull(reported_trait) %>%
  unique()

total_lymphocyte_count_assocs <- get_gws_associations_for_lymphocyte_trait(total_lymphocyte_count_traits)


b_and_t_cell_assocs <- fuzzy_join_assocs(
  gws_b_cell_assocs,
  gws_t_cell_assocs,
  "b_cell",
  "t_cell"
)

total_lymphocyte_count_and_t_cell_assocs <- fuzzy_join_assocs(
  total_lymphocyte_count_assocs,
  gws_t_cell_assocs,
  "l_counts",
  "t_cell"
)

total_lymphocyte_count_and_b_cell_assocs <- fuzzy_join_assocs(
  total_lymphocyte_count_assocs,
  gws_b_cell_assocs,
  "l_counts",
  "b_cell"
)

write.table(b_and_t_cell_assocs, sep = "\t", file = snakemake@input$b_and_t_cell_assocs)
write.table(gws_t_cell_assocs, sep = "\t", file = snakemake@input$t_cell_assocs)
write.table(gws_b_cell_assocs, sep = "\t", file = snakemake@input$b_cell_assocs)
write.table(total_lymphocyte_count_assocs, sep = "\t", file = snakemake@input$lymphocyte_count_assocs)
write.table(total_lymphocyte_count_and_t_cell_assocs, sep = "\t", file = snakemake@input$lymphocyte_count_and_t_cell_assocs)
write.table(total_lymphocyte_count_and_b_cell_assocs, sep = "\t", file = snakemake@input$lymphocyte_count_and_b_cell_assocs)
