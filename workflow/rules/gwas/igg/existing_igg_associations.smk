rule collate_existing_igg_associations:
    input:
        ebi = "resources/gwas/igg/gwas-association-downloaded_2025-02-04-EFO_0004565.tsv",
#        epic = "results/harmonised_gwas/epic-igg/1000kb_gws_annotated_lead_snps.tsv",
        pietzner = "results/harmonised_gwas/pietzner-igg/1000kb_gws_annotated_lead_snps.tsv",
        gudjonsson = "results/harmonised_gwas/gudjonsson-igg/1000kb_gws_annotated_lead_snps.tsv",
        eldjarn = "results/harmonised_gwas/eldjarn-igg/1000kb_gws_annotated_lead_snps.tsv",
        scepanovic = "results/harmonised_gwas/scepanovic-igg/1000kb_gws_annotated_lead_snps.tsv",
        dennis = "results/harmonised_gwas/dennis-igg/1000kb_gws_annotated_lead_snps.tsv",
    output:
        "results/igg_meta/existing_associations.tsv"
    localrule: True
    script: script_path("gwas/igg_meta/collate_existing_associations.R")

use rule merge_iga_meta_lead_snps_with_existing_associations as merge_igg_meta_lead_snps_with_existing_associations with:
    input:
        existing = "results/igg_meta/existing_associations.tsv",
        meta = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv"
