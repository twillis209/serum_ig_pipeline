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

rule draw_stacked_ig_loci_manhattans:
    input:
        sumstats = rules.combine_ighkl_regions_from_gwas.output,
        edb = rules.download_ensembl_db.output
    output:
        "results/gwas/ighkl/{isotype}_{ighkl_locus}.{ext}"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_stacked_ighkl_plots.R")

rule draw_ighkl_manhattans:
    input:
        expand(rules.draw_stacked_ig_loci_manhattans.output[0], isotype = ["iga", "igm", "igg"], ighkl_locus = ["igh", "igk", "igl"], ext = ["pdf"])

rule draw_igh_iga_associations:
    input:
        sumstats = rules.combine_ighkl_regions_from_gwas.output,
        edb = rules.download_ensembl_db.output
    output:
        "results/paper/figures/igh_iga_associations.png"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_igh_iga_associations.R")

rule draw_igh_igg_associations:
    input:
        sumstats = rules.combine_ighkl_regions_from_gwas.output,
        edb = rules.download_ensembl_db.output
    output:
        "results/paper/figures/igh_igg_associations.png"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_igh_igg_associations.R")
