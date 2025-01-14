use rule merge_iga_gwas as merge_iga_and_igg_meta with:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/with_gudjonsson/with_eldjarn/meta.tsv.gz",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/with_gudjonsson/with_eldjarn/meta.tsv.gz"
    output:
        "results/merged_gwas/iga-meta_and_igg-meta/{join,inner}/{variant_set}/merged.tsv.gz"