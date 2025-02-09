rule igh_associations_table:
    input:
        "results/ig/1000kb_study_and_meta_gws_igh_lead_snps.tsv"
    output:
        "results/paper/tables/igh_associations.tsv"
    localrule: True
    script: script_path("paper/tables/igh_associations.R")

rule h2_and_rg_estimates:
    input:
        "results/ig/h2_estimates.tsv",
        "results/ig/rg_estimates.tsv"

