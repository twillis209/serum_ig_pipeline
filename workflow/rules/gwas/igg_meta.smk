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
