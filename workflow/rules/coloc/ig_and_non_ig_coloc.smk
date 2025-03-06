checkpoint merge_ig_and_non_ig_lead_snps:
    input:
        ig = branch(evaluate("{isotype}"),
                    cases = {
                        "iga" : "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igg" : "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
                        "igm" : "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
                    }
        ),
        non_ig = "results/harmonised_gwas/{non_ig}/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/coloc/{isotype}_and_{non_ig}/merged_lead_snps.tsv",
    params:
        window = 1e6
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("coloc/merge_ig_and_non_ig_lead_snps.R")

rule subset_ig_and_non_ig_pair_for_coloc:
    input:
        lead_snps = "results/coloc/{isotype}_and_{non_ig}/merged_lead_snps.tsv",
        merged = "results/merged_gwas/{isotype}-meta_and_{non_ig}/inner/with_mhc/merged.tsv.gz"
    output:
        "results/coloc/{isotype}_and_{non_ig}/{first_rsid}_and_{second_rsid}/sumstats.tsv"
    params:
        flank = 250000
    threads: 8
    conda: env_path("global.yaml")
    script: script_path("coloc/subset_ig_and_non_ig_pair_for_locus.R")
