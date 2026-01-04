rule combine_ighkl_regions_from_gwas:
    input:
        iga = rules.merge_iga_gwas.output,
        iga_meta = rules.run_iga_meta_analysis.output,
        igg = rules.merge_igg_gwas.output,
        igg_meta = rules.run_igg_meta_analysis.output,
        igm = rules.merge_igm_gwas.output,
        igm_meta = rules.run_igm_meta_analysis.output,
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
