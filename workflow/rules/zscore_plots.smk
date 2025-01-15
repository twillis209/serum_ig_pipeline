rule draw_zscore_plot_for_ig_gwas:
    input:
        "results/{isotype}_meta/merged.tsv.gz"
    output:
        "results/{isotype}_meta/zscore_plots/{study_a}_{study_b}.png"
    params:
        beta_a = lambda w: f"{config.get('beta_col')}.{w.study_a}",
        beta_b = lambda w: f"{config.get('beta_col')}.{w.study_b}",
        se_a = lambda w: f"{config.get('se_col')}.{w.study_a}",
        se_b = lambda w: f"{config.get('se_col')}.{w.study_b}"
    threads: 8
    resources:
        runtime = 20
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/draw_zscore_plot.R")

rule zscore_plots:
    input:
        [[f"results/{x}_meta/zscore_plots/{study_a}_{study_b}.png" for study_a, study_b in combinations(config.get(f'{x}_studies'), 2)] for x in ['iga', 'igm', 'igg']]
