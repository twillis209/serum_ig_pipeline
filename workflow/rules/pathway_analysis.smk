rule extract_gene_sets:
    input:
        iga = "results/paper/tables/iga_lead_snps.tsv",
        igg = "results/paper/tables/igg_lead_snps.tsv",
        igm = "results/paper/tables/igm_lead_snps.tsv",
        ig_coloc = "results/paper/tables/ig_coloc.tsv"
    output:
        "results/pathway_analysis/ig_lead_snp_genes.txt"
    conda: env_path("pathway.yml")
    script: script_path("pathway/prepare_inputs.R")

rule draw_pathway_heatmap:
    input:
        "results/pathway_analysis/ig_lead_snp_genes.txt"
    output:
        png = "results/pathway_analysis/combined_pathways_heatmap.png",
        tiff = "results/pathway_analysis/combined_pathways_heatmap.tiff"
    conda: env_path("pathway.yml")
    script: script_path("pathway/pathway_heatmap_plot.R")

rule draw_coloc_heatmap:
    input:
        "results/paper/tables/ig_and_non_ig_coloc.tsv"
    output:
        "results/pathway_analysis/ig_coloc_heatmap.png"
    conda: env_path("pathway.yml")
    script: script_path("pathway/coloc_heatmap_plot.R")

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
