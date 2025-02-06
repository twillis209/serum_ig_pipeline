use rule collate_existing_iga_associations as collate_existing_igm_associations with:
    input:
        ebi = "resources/gwas/igm/gwas-association-downloaded_2025-02-04-EFO_0004993.tsv",
        pietzner = "results/harmonised_gwas/pietzner-igm/1000kb_gws_annotated_lead_snps.tsv",
        gudjonsson = "results/harmonised_gwas/gudjonsson-igm/1000kb_gws_annotated_lead_snps.tsv",
        eldjarn = "results/harmonised_gwas/eldjarn-igm/1000kb_gws_annotated_lead_snps.tsv",
        scepanovic = "results/harmonised_gwas/scepanovic-igm/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/igm_meta/existing_associations.tsv"

use rule merge_iga_meta_lead_snps_with_existing_associations as merge_igm_meta_lead_snps_with_existing_associations with:
    input:
        existing = "results/igm_meta/existing_associations.tsv",
        meta = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv"
