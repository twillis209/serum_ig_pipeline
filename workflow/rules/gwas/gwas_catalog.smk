rule fetch_lymphocyte_associations:
    output:
        b_cell_assocs = "results/gwas_catalog/b_cell_associations.tsv",
        t_cell_assocs = "results/gwas_catalog/t_cell_associations.tsv",
        lymphocyte_count_assocs = "results/gwas_catalog/lymphocyte_count_associations.tsv",
        b_and_t_cell_assocs = "results/gwas_catalog/b_and_t_cell_associations.tsv",
        lymphocyte_count_and_t_cell_assocs = "results/gwas_catalog/lymphocyte_count_and_t_cell_associations.tsv",
        lymphocyte_count_and_b_cell_assocs = "results/gwas_catalog/lymphocyte_count_and_b_cell_associations.tsv",
    localrule: True
    conda: env_path("gwasrapidd.yaml")
    script: script_path("gwas/gwas_catalog/fetch_lymphocyte_associations.R")
