rule combine_ighkl_regions_from_gwas:
    input:
        iga = "results/iga_meta/merged.tsv.gz",
        iga_meta = "resources/harmonised_gwas/iga-meta.tsv.gz",
        igg = "results/igg_meta/merged.tsv.gz",
        igg_meta = "resources/harmonised_gwas/igg-meta.tsv.gz",
        igm = "results/igm_meta/merged.tsv.gz",
        igm_meta = "resources/harmonised_gwas/igm-meta.tsv.gz"
    output:
        "results/gwas/ighkl/combined_ighkl_regions.tsv.gz"
    params:
        flank = 1000000
    threads: 16
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/merge_ighkl_regions.R")
