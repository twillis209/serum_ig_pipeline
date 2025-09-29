# Figure 1
rule draw_stacked_ig_manhattans:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        png = "results/paper/figures/ig_manhattans.png",
        tiff = "results/paper/figures/ig_manhattans.tiff"
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

# Figure 2
rule draw_pathway_plots:
    input:
        "results/pathway_analysis/ig_lead_snp_genes.txt"
    output:
        combined_pathways_heatmap_png = "results/pub/figures/combined_pathways_heatmap.png",
        combined_pathways_enrichment_png = "results/pub/figures/combined_pathways_enrichment.png",
        combined_pathways_heatmap_tiff = "results/pub/figures/combined_pathways_heatmap.tiff",
        combined_pathways_enrichment_tiff = "results/pub/figures/combined_pathways_enrichment.tiff",
    localrule: True
    conda: env_path("pathway.yml")
    script: script_path("pathway/draw_pathway_plots.R")

# Supplementary figures
rule supp_figures:
    input:
        expand("results/coloc/iga_and_igm/{rsid_pair}/trimmed/lz_plots.tiff", rsid_pair = ["rs188468174_and_rs188468174", "rs6570530_and_rs7758383", "rs35969813_and_rs35629860", "rs2872516_and_rs4795397", "rs6885567_and_rs6879652", "rs10517538_and_rs1397934"]),
        "results/coloc/iga_and_igg/rs3803800_and_rs758641530/trimmed/lz_plots.tiff",
        "results/coloc/iga_and_igg/rs3803800_and_rs758641530/untrimmed/lz_plots.tiff",
        "results/coloc/raw_vs_trimmed_coloc_pp_h4_ig_pairs.tiff",
        "results/pub/figures/combined_pathways_heatmap.tiff"
