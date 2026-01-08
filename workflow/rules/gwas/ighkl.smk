rule combine_ighkl_regions_from_gwas:
    input:
        iga = "results/iga_meta/merged.tsv.gz",
        iga_meta = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz",
        igg = "results/igg_meta/merged.tsv.gz",
        igg_meta = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz",
        igm = "results/igm_meta/merged.tsv.gz",
        igm_meta = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz"
    output:
        "results/gwas/ighkl/{ighkl_flank}/combined_ighkl_regions.tsv.gz"
    params:
        flank = lambda w: int(w.ighkl_flank)
    threads: 16
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/merge_ighkl_regions.R")

rule draw_stacked_ig_loci_manhattans:
    input:
        sumstats = "results/gwas/ighkl/1000000/combined_ighkl_regions.tsv.gz",
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
        sumstats = "results/gwas/ighkl/1000000/combined_ighkl_regions.tsv.gz",
        edb = rules.download_ensembl_db.output
    output:
        "results/paper/figures/igh_iga_associations.png"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_igh_iga_associations.R")

rule draw_igh_igg_associations:
    input:
        sumstats = "results/gwas/ighkl/1000000/combined_ighkl_regions.tsv.gz",
        edb = rules.download_ensembl_db.output
    output:
        "results/paper/figures/igh_igg_associations.png"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_igh_igg_associations.R")

rule draw_igh_igm_associations:
    input:
        sumstats = "results/gwas/ighkl/1000000/combined_ighkl_regions.tsv.gz",
        edb = rules.download_ensembl_db.output
    output:
        "results/paper/figures/igh_igm_associations.png"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_igh_igm_associations.R")

rule draw_igk_igg_associations:
    input:
        sumstats = "results/gwas/ighkl/1000000/combined_ighkl_regions.tsv.gz",
        edb = rules.download_ensembl_db.output
    output:
        "results/paper/figures/igk_igg_associations.png"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/draw_igk_igg_associations.R")
