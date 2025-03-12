rule draw_locuszoomr_plot_for_interval:
    input:
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        "results/locuszoomr/{trait}/{tag}_chr{chrom}_{start}_{end}.png"
    threads: 8
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("gwas/locuszoomr/plot_locus.R")
