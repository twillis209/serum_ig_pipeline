rule combine_ighkl_regions_from_gwas:
    input:
        iga = "results/iga_meta/merged.tsv.gz",
        iga_meta = "results/iga_meta/meta.tsv.gz",
        igg = "results/igg_meta/merged.tsv.gz",
        igg_meta = "results/igg_meta/meta.tsv.gz",
        igm = "results/igm_meta/merged.tsv.gz",
        igm_meta = "results/igm_meta/meta.tsv.gz"
    output:
        "results/gwas/ighkl/combined_ighkl_regions.tsv.gz"
    params:
        flank = 1000000
    threads: 16
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/merge_ighkl_regions.R")
