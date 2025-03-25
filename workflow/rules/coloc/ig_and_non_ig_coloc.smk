def get_rsids_from_merged_lead_snps_for_ig_non_ig_pair(w):
    daf = pd.read_csv(checkpoints.merge_ig_and_non_ig_lead_snps.get(isotype = w.isotype, non_ig = w.non_ig).output[0], sep = '\t')

    return zip(daf[f"rsid.{w.isotype}"], daf[f"rsid.{w.non_ig}"])

checkpoint merge_ig_and_non_ig_lead_snps:
    input:
        ig = branch(evaluate("{isotype}"),
                    cases = {
                        "iga" : "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igg" : "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igm" : "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
                    }
        ),
        non_ig = "results/harmonised_gwas/{non_ig}/{non_ig_window}_gws_annotated_lead_snps.tsv"
    output:
        "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/merged_lead_snps.tsv"
    params:
        window = 1e6,
        loci_to_drop = ['igh', 'igk', 'igl']
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/merge_ig_and_non_ig_lead_snps.R")

rule subset_ig_and_non_ig_pair_for_coloc:
    input:
        lead_snps = "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/merged_lead_snps.tsv",
        merged = "results/merged_gwas/{isotype}-meta_and_{non_ig}/inner/with_mhc/merged.tsv.gz"
    output:
        "results/coloc/{isotype}_and_{non_ig_window}_{non_ig}/{non_ig_window}/{isotype_rsid}_and_{non_ig_rsid}/sumstats.tsv"
    params:
        flank = 250000
    threads: 8
    conda: env_path("global.yaml")
    script: script_path("coloc/subset_ig_and_non_ig_pair_for_locus.R")

rule run_coloc_for_ig_and_non_ig_pair:
    input:
        "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/{isotype_rsid}_and_{non_ig_rsid}/sumstats.tsv"
    output:
        rds = "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/{isotype_rsid}_and_{non_ig_rsid}/coloc.rds",
        tsv = "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/{isotype_rsid}_and_{non_ig_rsid}/coloc.tsv"
    conda: env_path("coloc.yaml")
    script: script_path("coloc/run_coloc_for_ig_and_non_ig_pair.R")

use rule run_coloc_for_all_ig_pairs as run_coloc_for_all_snps_for_ig_and_non_ig_pair with:
    input:
        lambda w: [f"results/coloc/{{isotype}}_and_{{non_ig}}/{{non_ig_window}}/{isotype_rsid}_and_{non_ig_rsid}/coloc.tsv" for isotype_rsid, non_ig_rsid in get_rsids_from_merged_lead_snps_for_ig_non_ig_pair(w)]
    output:
        "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/results.tsv"

rule draw_locuszoom_plots_for_all_snps_for_ig_and_non_ig_pair:
    input:
        lambda w: [f"results/coloc/{{isotype}}_and_{{non_ig}}/{{non_ig_window}}/{isotype_rsid}_and_{non_ig_rsid}/lz_plots.png" for isotype_rsid, non_ig_rsid in get_rsids_from_merged_lead_snps_for_ig_non_ig_pair(w)]
    output:
        "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/lz_plots.done"
    shell: "touch {output}"

rule add_genes_to_ig_and_non_ig_coloc_pair:
    input:
        ig = branch(evaluate("{isotype}"),
                    cases = {
                        "iga" : "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igg" : "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igm" : "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
                    }
                    ),
        # I think 3Mb windows should be ok for asthma and lymphocyte-counts
        non_ig = "results/harmonised_gwas/{non_ig}/{non_ig_window}_gws_annotated_lead_snps.tsv",
        coloc = "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/results.tsv"
    output:
        "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/results_with_genes.tsv"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/add_genes_to_ig_and_non_ig_coloc_pair.R")

use rule draw_locuszoomr_plot_for_coloc_ig_pair as draw_locuszoomr_plot_for_coloc_ig_and_non_ig_pair with:
    input:
        sumstats = "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/{first_rsid}_and_{second_rsid}/sumstats.tsv",
        coloc = "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/{first_rsid}_and_{second_rsid}/coloc.tsv"
    output:
        "results/coloc/{isotype}_and_{non_ig}/{non_ig_window}/{first_rsid}_and_{second_rsid}/lz_plots.png"
