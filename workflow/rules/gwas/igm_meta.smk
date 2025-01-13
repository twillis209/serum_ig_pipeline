use rule merge_iga_gwas as merge_igm_gwas with:
    input:
        epic = "results/restandardised_gwas/epic-igm.tsv.gz",
        scepanovic = "results/restandardised_gwas/scepanovic-igm.tsv.gz",
        pietzner = "resources/harmonised_gwas/pietzner-igm.tsv.gz",
        gudjonsson = "resources/harmonised_gwas/gudjonsson-igm.tsv.gz",
        eldjarn = "results/restandardised_gwas/eldjarn-igm.tsv.gz",
    output:
        "results/igm_meta/merged.tsv.gz"

use rule run_iga_meta_analysis as run_igm_meta_meta_analysis with:
    input:
        "results/igm_meta/merged.tsv.gz"
    output:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"

checkpoint distance_clump_igm_meta:
    input:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    output:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = lambda w: False,
        index_threshold = lambda w: 5e-8 if w.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
        snps_to_ignore = []
    threads: 16
    resources:
        runtime = 5
    conda: env_path("global.yaml")
    script: script_path("gwas/distance_clump.R")

use rule annotate_lead_snps as annotate_igm_meta_lead_snps with:
    input:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_lead_snps.tsv"
    output:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_annotated_lead_snps.tsv"

use rule draw_manhattan_with_lead_snp_annotation as draw_igm_meta_manhattan with:
    input:
        gwas = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz",
        lead_snps = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/manhattan.png"
    params:
        title = '',
        width = 6,
        height = 4

use rule make_plink_range as make_plink_range_for_igm_meta with:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    output:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt"

use rule subset_reference as subset_reference_for_igm_meta with:
    input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        temp(multiext("results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"))

use rule calculate_human_default_taggings as calculate_human_default_taggings_for_igm_meta with:
    input:
        multiext("results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output:
        tagging_file = temp("results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging")

use rule process_sum_stats as process_sum_stats_for_igm_meta with:
    input:
        gwas_file = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz",
        range_file = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt",
    output:
        "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/procd.assoc"

use rule estimate_h2_with_human_default as estimate_h2_with_human_default_for_igm_meta with:
    input:
        gwas = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/procd.assoc",
        tagging_file = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging"
    output:
        multiext("results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
