rule merge_iga_gwas:
    input:
        epic = "results/restandardised_gwas/epic-iga.tsv.gz",
        liu = "results/restandardised_gwas/liu-iga.tsv.gz",
        dennis = "results/restandardised_gwas/dennis-iga.tsv.gz",
        pietzner = "resources/harmonised_gwas/pietzner-iga.tsv.gz",
        eldjarn = "results/restandardised_gwas/eldjarn-iga.tsv.gz",
        gudjonsson = "resources/harmonised_gwas/gudjonsson-iga.tsv.gz",
        scepanovic = "results/restandardised_gwas/scepanovic-iga.tsv.gz"
    output:
        "results/iga_meta/merged.tsv.gz"
    threads: 24
    resources:
        runtime = 120
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/merge_named_gwas_inputs.R")

rule run_iga_meta_analysis:
    input:
        "results/iga_meta/merged.tsv.gz"
    output:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    threads: 16
    resources:
        runtime = 30
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/run_meta_analysis.R")

rule drop_selected_loci_from_iga_meta_analysis:
    input:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    output:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/filtered_meta.tsv.gz"
    threads: 16
    resources:
        runtime = 30
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/drop_loci_from_meta_analysis.R")

checkpoint distance_clump_iga_meta:
    input:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz"
    output:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_lead_snps.tsv"
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

use rule annotate_lead_snps as annotate_iga_meta_lead_snps with:
    input:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_lead_snps.tsv"
    output:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{window_size}_{threshold}_annotated_lead_snps.tsv"

use rule draw_manhattan_with_lead_snp_annotation as draw_iga_meta_manhattan with:
    input:
        gwas = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/meta.tsv.gz",
        lead_snps = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/2000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/manhattan.png"
    params:
        title = '',
        width = 6,
        height = 4,
        ylim = [0, 1e-80],

use rule make_plink_range as make_plink_range_for_iga_meta with:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/filtered_meta.tsv.gz"
    output:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt"

use rule subset_reference as subset_reference_for_iga_meta with:
    input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        temp(multiext("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"))

use rule calculate_human_default_taggings as calculate_human_default_taggings_for_iga_meta with:
    input:
        multiext("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output:
        tagging_file = temp("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging")
    log:
        log_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging.log"

use rule process_sum_stats as process_sum_stats_for_iga_meta with:
    input:
        gwas_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/filtered_meta.tsv.gz",
        range_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/matching_ids.txt",
    output:
        temp("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/procd.assoc")
    params:
        N = lambda w: get_combined_sample_size_for_ig_meta(w, 'iga')

use rule estimate_h2_with_human_default as estimate_h2_with_human_default_for_iga_meta with:
    input:
        gwas = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/procd.assoc",
        tagging_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/merged.tagging"
    output:
        multiext("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/sumher.log"

rule draw_igh_locus_for_iga_datasets:
    input:
        "results/iga_meta/merged.tsv.gz"
    output:
        "results/iga_meta/locuszoomr/igh.pdf"
    params:
        p_value_cols = [f"p_value.{x}" for x in config['iga_studies']],
        chrom = config['loci']['igh']['chrom'],
        start_pos = config['loci']['igh']['start'] - 250000,
        stop_pos = config['loci']['igh']['stop'] + 250000
    threads: 12
    conda: env_path("global.yaml")
    script: script_path("gwas/locuszoomr/plot_locus.R")

#
#use rule subset_summary_statistics_about_variant as subset_summary_statistics_about_variant_for_iga_meta with:
#    input:
#        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/meta_prescreen.tsv.gz"
#    output:
#        sum_stats = temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz"),
#        ids = temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ids.txt")
#
#use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_iga_meta_ld_matrix with:
#    input:
#        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
#        ids = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ids.txt"
#    output:
#        temp(multiext("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam"))
#    params:
#        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
#        out_stem = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset"
#
#use rule calculate_ld_for_subset_about_variant as calculate_ld_for_subset_about_iga_meta_variant with:
#    input:
#        multiext("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam")
#    output:
#        temp(multiext("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset.phased", ".vcor2", ".vcor2.vars"))
#    params:
#        stem = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset",
#        variant_id = lambda w: w.variant_id.replace('_', ':')
#
#use rule fetch_rsids_for_sum_stats_about_igad_meta_variant as fetch_rsids_for_sum_stats_about_iga_meta_variant with:
#    input:
#        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz"
#    output:
#        temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats_with_rsids.tsv.gz")
#
#use rule merge_sumstats_with_r2 as merge_iga_meta_sumstats_with_r2 with:
#    input:
#        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz",
#        ld = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset.phased.vcor2",
#        ld_vars = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset.phased.vcor2.vars"
#    output:
#        temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ld_friends.tsv.gz")
#
#use rule draw_locuszoom_plot_without_r2 as draw_locuszoom_plot_without_r2_for_iga_meta with:
#    input:
#        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz",
#        rsIDs = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.rsIDs"
#    output:
#        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/locuszoom_sans_r2_{gene_track}.png"
#
#use rule draw_locuszoom_plot_with_r2 as draw_locuszoom_plot_with_r2_for_iga_meta with:
#    input:
#        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ld_friends.tsv.gz",
#        rsIDs = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.rsIDs"
#    output:
#        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/locuszoom_with_r2_{gene_track}.png"
#
#rule plot_all_meta_gws_loci:
#    input:
#        [f"results/iga_meta/with_decode/with_dennis/prescreen/gws/locus_plots/{locus}/index_snp.png" for locus in config.get('iga').get('all')]
#    output:
#        "results/iga_meta/with_decode/with_dennis/prescreen/gws/locus_plots/meta_gws.done"
#    localrule: True
#    shell: "touch {output}"
#
