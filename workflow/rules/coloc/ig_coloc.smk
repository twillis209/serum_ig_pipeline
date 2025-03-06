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
        window = 1e6
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
        rds = "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/coloc.rds",
        tsv = "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/coloc.tsv"
    params:
        first_isotype_max_n = lambda w: config.get('gwas_datasets').get(w.first_isotype),
        second_isotype_max_n = lambda w: config.get('gwas_datasets').get(w.second_isotype)
    conda: env_path("coloc.yaml")
    script: script_path("coloc/run_coloc.R")

rule draw_coloc_plots_for_ig_pair:
    input:
        "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/sumstats.tsv"
    output:
        "results/coloc/{first_isotype}_and_{second_isotype}/{first_rsid}_and_{second_rsid}/plots.png"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/draw_coloc_plots.R")

rule run_coloc_for_all_ig_pairs:
    input:
        iga_igg = lambda w: [f"results/coloc/iga_and_igg/{first_rsid}_and_{second_rsid}/coloc.tsv" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igg')],
        igg_igm = lambda w: [f"results/coloc/igg_and_igm/{first_rsid}_and_{second_rsid}/coloc.tsv" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'igg', isotype_b = 'igm')],
        iga_igm = lambda w: [f"results/coloc/iga_and_igm/{first_rsid}_and_{second_rsid}/coloc.tsv" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igm')]
    output:
        "results/coloc/all_ig_pairs.tsv"
    localrule: True
    run:
        dafs = [pd.read_csv(x, sep = '\t') for x in input]

        pd.concat(dafs).to_csv(output[0], sep = '\t', index = False)

rule add_gene_to_all_ig_coloc_pairs:
    input:
        all_pairs = "results/coloc/all_ig_pairs.tsv",
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/coloc/all_ig_pairs_with_genes.tsv"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/add_genes_to_ig_coloc_pairs.R")

rule draw_coloc_plots_for_all_ig_pairs:
    input:
        iga_igg = lambda w: [f"results/coloc/iga_and_igg/{first_rsid}_and_{second_rsid}/plots.png" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igg')],
        igg_igm = lambda w: [f"results/coloc/igg_and_igm/{first_rsid}_and_{second_rsid}/plots.png" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'igg', isotype_b = 'igm')],
        iga_igm = lambda w: [f"results/coloc/iga_and_igm/{first_rsid}_and_{second_rsid}/plots.png" for first_rsid, second_rsid in get_rsids_from_merged_lead_snps(w, isotype_a = 'iga', isotype_b = 'igm')]
