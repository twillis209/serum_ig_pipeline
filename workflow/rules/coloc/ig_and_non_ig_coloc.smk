def get_rsids_from_candidate_lead_snps_for_ig_non_ig_pair(w):
    daf = pd.read_csv(checkpoints.tabulate_ig_and_non_ig_coloc_candidates.get(isotype = w.isotype, non_ig = w.non_ig).output[0], sep = '\t')

    return daf[daf['is_candidate'] == True]["rsid"]

checkpoint tabulate_ig_and_non_ig_coloc_candidates:
    input:
        ig = branch(evaluate("{isotype}"),
                    cases = {
                        "iga" : "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igg" : "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igm" : "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
                    }
        ),
        merged = "results/merged_gwas/{isotype}-meta_and_{non_ig}/inner/with_mhc/merged.tsv.gz"
    output:
        "results/coloc/{isotype}_and_{non_ig}/coloc_candidate_lead_snps.tsv"
    params:
        window = 1e6
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/tabulate_ig_and_non_ig_coloc_candidates.R")

rule subset_ig_and_non_ig_pair_for_coloc:
    input:
        lead_snps = "results/coloc/{isotype}_and_{non_ig}/coloc_candidate_lead_snps.tsv",
        merged = "results/merged_gwas/{isotype}-meta_and_{non_ig}/inner/with_mhc/merged.tsv.gz"
    output:
        "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/sumstats.tsv"
    params:
        flank = 500000
    threads: 8
    resources:
        runtime = 10
    group: "coloc"
    conda: env_path("global.yaml")
    script: script_path("coloc/subset_ig_and_non_ig_pair_for_locus.R")

rule run_coloc_for_ig_and_non_ig_pair:
    input:
        "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/sumstats.tsv"
    output:
        rds = "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/coloc.rds",
        tsv = "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/coloc.tsv"
    resources:
        runtime = 10
    group: "coloc"
    conda: env_path("coloc.yaml")
    script: script_path("coloc/run_coloc_for_ig_and_non_ig_pair.R")

rule add_genes_and_r2_to_ig_and_non_ig_coloc_pair:
    input:
        ig = "results/coloc/{isotype}_and_{non_ig}/coloc_candidate_lead_snps.tsv",
        coloc = "results/coloc/{isotype}_and_{non_ig}/results.tsv"
    output:
        "results/coloc/{isotype}_and_{non_ig}/results_with_genes_and_r2.tsv"
    resources:
        ldlink_calls = 1
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/add_genes_and_r2_to_ig_and_non_ig_coloc_pair.R")

use rule run_coloc_for_all_ig_pairs as run_coloc_for_all_snps_for_ig_and_non_ig_pair with:
    input:
        lambda w: [f"results/coloc/{{isotype}}_and_{{non_ig}}/{isotype_rsid}/coloc.tsv" for isotype_rsid in get_rsids_from_candidate_lead_snps_for_ig_non_ig_pair(w)]
    output:
        "results/coloc/{isotype}_and_{non_ig}/results.tsv"
    params:
        cols = ["nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "first_trait", "second_trait", "ig_snp", "non_ig_snp", "min_p.first", "min_p.second", "pearson.cor", "ig_snp_effect_ratio", "non_ig_snp_effect_ratio"]

rule draw_locuszoom_plots_for_all_snps_for_ig_and_non_ig_pair:
    input:
        lambda w: [f"results/coloc/{{isotype}}_and_{{non_ig}}/{isotype_rsid}/lz_plots.png" for isotype_rsid in get_rsids_from_candidate_lead_snps_for_ig_non_ig_pair(w)] or []
    output:
        "results/coloc/{isotype}_and_{non_ig}/lz_plots.done"
    localrule: True
    shell: "touch {output}"

use rule draw_locuszoomr_plot_for_coloc_ig_pair as draw_locuszoomr_plot_for_coloc_ig_and_non_ig_pair with:
    input:
        sumstats = "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/sumstats.tsv",
        coloc = "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/coloc.tsv"
    output:
        "results/coloc/{isotype}_and_{non_ig}/{isotype_rsid}/lz_plots.png"
    params:

rule ig_and_non_ig_coloc_tables_and_locuszoomr_plots:
    input:
        expand("results/coloc/{isotype}_and_{non_ig}/{filetype}",
               isotype = ["igg", "iga", "igm"],
               non_ig = config['imds'],
               filetype = ["lz_plots.done", "results_with_genes_and_r2.tsv"])
