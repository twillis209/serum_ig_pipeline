rule igh_associations_table:
    input:
        epic_iga = "results/harmonised_gwas/epic-iga/2000kb_gws_igh_lead_snps.tsv",
        epic_igg = "results/harmonised_gwas/epic-igg/2000kb_gws_igh_lead_snps.tsv",
        pietzner_iga = "results/harmonised_gwas/pietzner-iga/2000kb_gws_igh_lead_snps.tsv",
        pietzner_igg = "results/harmonised_gwas/pietzner-igg/2000kb_gws_igh_lead_snps.tsv",
        eldjarn_iga = "results/harmonised_gwas/eldjarn-iga/2000kb_gws_igh_lead_snps.tsv",
        eldjarn_igg = "results/harmonised_gwas/eldjarn-igg/2000kb_gws_igh_lead_snps.tsv",
        eldjarn_igm = "results/harmonised_gwas/eldjarn-igm/2000kb_gws_igh_lead_snps.tsv"
    output:
        "results/paper/tables/igh_associations.tsv"
    localrule: True
    script: script_path("paper/tables/igh_associations.R")
