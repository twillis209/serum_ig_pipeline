rule extract_gene_sets:
    input:
        iga = "results/paper/tables/iga_lead_snps.tsv",
        igg = "results/paper/tables/igg_lead_snps.tsv",
        igm = "results/paper/tables/igm_lead_snps.tsv",
        ig_coloc = "results/paper/tables/ig_coloc.tsv"
    output:
        "results/pathway_analysis/ig_lead_snp_genes.txt"

rule draw_pathway_plots:
    input:
        "results/pathway_analysis/ig_lead_snp_genes.txt"
    output:
        combined_pathways_heatmap = "results/pub/figures/combined_pathways_heatmap.png",
        combined_pathways_enrichment_plot = "results/pub/figures/combined_pathways_heatmap.png",
    localrule: True
    conda: env_path("pathway.yml")
    script: script_path("pathway_analysis/draw_pathway_plots.R")
