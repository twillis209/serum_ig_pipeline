rule run_mr_on_ig_and_non_ig:
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
        tsv = "results/mr/{isotype}_and_{non_ig}/mr.tsv",
        png = "results/mr/{isotype}_and_{non_ig}/mr.png",
        png_ivw = "results/mr/{isotype}_and_{non_ig}/mr_ivw.png"
    params:
        snps_to_exclude = []
    threads: 8
    localrule: True
    conda: env_path("mr.yaml")
    script: script_path("mr/run_mr.R")
