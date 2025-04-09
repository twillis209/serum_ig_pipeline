def get_rsids_from_merged_lead_snps(w, isotype_a, isotype_b):
    daf = pd.read_csv(checkpoints.merge_ig_lead_snps.get().output[f"{isotype_a}_{isotype_b}"], sep = '\t')

    return zip(daf[f"rsid.{isotype_a}"], daf[f"rsid.{isotype_b}"])

checkpoint merge_ig_lead_snps:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
    output:
        iga_igg = "results/coloc/iga_and_igg/merged_lead_snps.tsv",
        iga_igm = "results/coloc/iga_and_igm/merged_lead_snps.tsv",
        igg_igm = "results/coloc/igg_and_igm/merged_lead_snps.tsv"
    params:
        window = 1e6,
        loci_to_drop = ['igh', 'igk', 'igl']
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/merge_ig_lead_snps.R")

rule subset_ig_pair_for_coloc:
    input:
        lead_snps = "results/coloc/{first_isotype}_and_{second_isotype}/merged_lead_snps.tsv",
        merged = "results/merged_gwas/{first_isotype}_and_{second_isotype}/inner/with_mhc/merged.tsv.gz"
    output:
        "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/sumstats.tsv"
    params:
        flank = 250000
    threads: 8
    conda: env_path("global.yaml")
    script: script_path("coloc/subset_ig_pair_for_locus.R")

rule run_coloc_for_ig_pair:
    input:
        "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/sumstats.tsv"
    output:
        rds = "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/{trim}/coloc.rds",
        tsv = "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/{trim}/coloc.tsv"
    params:
        first_isotype_max_n = lambda w: config.get('gwas_datasets').get(w.first_isotype),
        second_isotype_max_n = lambda w: config.get('gwas_datasets').get(w.second_isotype)
    conda: env_path("coloc.yaml")
    script: script_path("coloc/run_coloc.R")

rule draw_locuszoomr_plot_for_coloc_ig_pair:
    input:
        sumstats = "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/sumstats.tsv",
        coloc = "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/{trim}/coloc.tsv"
    output:
        "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/{trim}/lz_plots.png"
    params:
        first_isotype_max_n = lambda w: config.get('gwas_datasets').get(w.first_isotype),
        second_isotype_max_n = lambda w: config.get('gwas_datasets').get(w.second_isotype)
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/draw_locuszoomr_plots.R")

rule run_coloc_for_all_ig_pairs:
    input:
        iga_igg = lambda w: [f"results/coloc/iga_and_igg/{first_rsid}_and_{second_rsid}/{trim}/coloc.tsv" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igg') for trim in ['trimmed', 'untrimmed']],
        igg_igm = lambda w: [f"results/coloc/igg_and_igm/{first_rsid}_and_{second_rsid}/{trim}/coloc.tsv" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'igg', isotype_b = 'igm') for trim in ['trimmed', 'untrimmed']],
        iga_igm = lambda w: [f"results/coloc/iga_and_igm/{first_rsid}_and_{second_rsid}/{trim}/coloc.tsv" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igm') for trim in ['trimmed', 'untrimmed']]
    output:
        "results/coloc/all_ig_pairs.tsv"
    localrule: True
    run:
        dafs = [pd.read_csv(x, sep = '\t') for x in input]

        pd.concat(dafs).to_csv(output[0], sep = '\t', index = False)

rule add_gene_and_r2_to_all_ig_coloc_pairs:
    input:
        all_pairs = "results/coloc/all_ig_pairs.tsv",
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/coloc/all_ig_pairs_with_genes_and_r2.tsv"
    resources:
        ldlink_calls = 1
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/add_genes_and_r2_to_ig_coloc_pairs.R")

rule draw_locuszoomr_plots_for_all_ig_pairs:
    input:
        iga_igg = lambda w: [f"results/coloc/iga_and_igg/{first_rsid}_and_{second_rsid}/{trim}/lz_plots.png" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igg') for trim in ['trimmed', 'untrimmed']],
        igg_igm = lambda w: [f"results/coloc/igg_and_igm/{first_rsid}_and_{second_rsid}/{trim}/lz_plots.png" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'igg', isotype_b = 'igm') for trim in ['trimmed', 'untrimmed']],
        iga_igm = lambda w: [f"results/coloc/iga_and_igm/{first_rsid}_and_{second_rsid}/{trim}/lz_plots.png" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igm') for trim in ['trimmed', 'untrimmed']]

rule draw_raw_vs_trimmed_coloc_pp_for_ig_pairs:
    input:
        "results/coloc/all_ig_pairs_with_genes_and_r2.tsv"
    output:
        h4 = "results/coloc/raw_vs_trimmed_coloc_pp_h4_ig_pairs.png"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/plot_raw_vs_trimmed_coloc_pp.R")
