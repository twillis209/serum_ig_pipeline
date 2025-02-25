rule fetch_coordinates_for_iuis_list:
    input:
        "resources/iei/IUIS-IEI-list-for-web-site-July-2024V2.tsv"
    output:
        sanitised = "results/iei/sanitised_iei_table.tsv"
    localrule: True
    conda: env_path("biomart.yaml")
    script: script_path("iei/add_gene_coordinates.R")
