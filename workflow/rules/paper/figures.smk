rule draw_stacked_ig_manhattans:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "results/paper/figures/ig_manhattans.png"
    params:
        width = 12,
        height = 14,
        iga_ylim = [1, 1e-60],
        igg_ylim = [1, 1e-40],
        igm_ylim = [1, 1e-60]
    threads: 12
    resources:
        runtime = 20
    conda: env_path("global.yaml")
    script: script_path("paper/figures/draw_stacked_ig_manhattans.R")

rule ig_ig_coloc_plots:
    input:
        "results/coloc/iga_and_igm/rs188468174_and_rs188468174/trimmed/lz_plots.png",
        "results/coloc/iga_and_igm/rs6570530_and_rs7758383/trimmed/lz_plots.png",
        "results/coloc/iga_and_igm/rs35969813_and_rs35629860/trimmed/lz_plots.png",
        "results/coloc/iga_and_igm/rs2872516_and_rs4795397/trimmed/lz_plots.png",
        "results/coloc/iga_and_igm/rs6885567_and_rs6879652/trimmed/lz_plots.png",
        "results/coloc/iga_and_igm/rs10517538_and_rs1397934/trimmed/lz_plots.png",
        "results/coloc/iga_and_igg/rs3803800_and_rs758641530/trimmed/lz_plots.png",
        "results/coloc/iga_and_igg/rs3803800_and_rs758641530/untrimmed/lz_plots.png",
        "results/coloc/raw_vs_trimmed_coloc_pp_h4_ig_pairs.png"
    output:
        "export/coloc/iga_and_igm/rs188468174_and_rs188468174/trimmed/lz_plots.png",
        "export/coloc/iga_and_igm/rs6570530_and_rs7758383/trimmed/lz_plots.png",
        "export/coloc/iga_and_igm/rs35969813_and_rs35629860/trimmed/lz_plots.png",
        "export/coloc/iga_and_igm/rs2872516_and_rs4795397/trimmed/lz_plots.png",
        "export/coloc/iga_and_igm/rs6885567_and_rs6879652/trimmed/lz_plots.png",
        "export/coloc/iga_and_igm/rs10517538_and_rs1397934/trimmed/lz_plots.png",
        "export/coloc/iga_and_igg/rs3803800_and_rs758641530/trimmed/lz_plots.png",
        "export/coloc/iga_and_igg/rs3803800_and_rs758641530/untrimmed/lz_plots.png",
        "export/coloc/raw_vs_trimmed_coloc_pp_h4_ig_pairs.png"
    localrule: True
    run:
        for i,x in enumerate(input):
            y = output[i]
            shell(f"cp {x} {y}")


