rule fetch_b_cell_associations:
    output:
        "results/gwas_catalog/b_cell_associations.tsv"
    localrule: True
    conda: env_path("gwasrapidd.yaml"
    script: script_path("gwas/gwas_catalog/fetch_b_cell_associations.R")
