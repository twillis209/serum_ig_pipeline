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

rule iga_lead_snps:
    input:
        lead = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv",
        novel = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv",
    output:
        "results/paper/tables/iga_lead_snps.tsv"
    localrule: True
    script: script_path("paper/tables/ig_lead_snp_table.R")
