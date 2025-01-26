use rule merge_iga_gwas as merge_igg_gwas with:
    input:
        epic = "results/restandardised_gwas/epic-igg.tsv.gz",
        dennis = "results/restandardised_gwas/dennis-igg.tsv.gz",
        scepanovic = "results/restandardised_gwas/scepanovic-igg.tsv.gz",
        pietzner = "resources/harmonised_gwas/pietzner-igg.tsv.gz",
        gudjonsson = "results/restandardised_gwas/gudjonsson-igg.tsv.gz",
        eldjarn = "resources/harmonised_gwas/eldjarn-igg.tsv.gz",
    output:
        "results/igg_meta/merged.tsv.gz"

use rule run_iga_meta_analysis as run_igg_meta_meta_analysis with:
    input:
        "results/igg_meta/merged.tsv.gz"
    output:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"

checkpoint distance_clump_igg_meta:
    input:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    output:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = lambda w: False,
        index_threshold = lambda w: 5e-8 if w.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
        snps_to_ignore = []
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/distance_clump.R")

use rule annotate_lead_snps as annotate_igg_meta_lead_snps with:
    input:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_lead_snps.tsv"
    output:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_annotated_lead_snps.tsv"

use rule draw_manhattan_with_lead_snp_annotation as draw_igg_meta_manhattan with:
    input:
        gwas = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz",
        lead_snps = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/manhattan.png"
    params:
        title = '',
        width = 6,
        height = 4

use rule make_plink_range as make_plink_range_for_igg_meta with:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    output:
        "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt"

use rule subset_reference as subset_reference_for_igg_meta with:
    input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        temp(multiext("results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"))

use rule calculate_human_default_taggings as calculate_human_default_taggings_for_igg_meta with:
    input:
        multiext("results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output:
        tagging_file = temp("results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging")
    log:
        log_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging.log"

use rule process_sum_stats as process_sum_stats_for_igg_meta with:
    input:
        gwas_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz",
        range_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt",
    output:
        temp("results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/procd.assoc")
    params:
        N = lambda w: get_combined_sample_size_for_ig_meta(w, 'igg')

use rule estimate_h2_with_human_default as estimate_h2_with_human_default_for_igg_meta with:
    input:
        gwas = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/procd.assoc",
        tagging_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging"
    output:
        multiext("results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = "results/igg_meta/{epic_inclusion}/{dennis_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/sumher.log"
