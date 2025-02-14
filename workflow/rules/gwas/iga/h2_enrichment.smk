rule write_out_ighkl_annotations_for_iga_meta:
    input:
        "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt"
    output:
        igh = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_1",
        igk = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_2",
        igl = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_3",
        base = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_4"
    threads: 8
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("ldsc_and_sumher/write_out_ighkl_annotations.R")

rule calculate_human_default_taggings_for_ighkl_annotations_and_iga_meta:
    input:
        multiext("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/merged", ".bed", ".bim", ".fam"),
        igh = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/ighkl_1",
        igk = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/ighkl_2",
        igl = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/ighkl_3",
        base = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/ighkl_4"
    output :
        tagging_file = temp("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_annotations/merged.tagging"),
    log:
        log_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/ighkl_annotations/merged.tagging.log"
    params:
        in_stem = subpath(input[0], strip_suffix = ".bed"),
        out_stem = subpath(output[0], strip_suffix = ".tagging"),
        partition_stem = subpath(input.igh, strip_suffix = "1")
    threads: 8
    resources:
        runtime = 30
    group: "sumher"
    shell:
        # NB: weightings now only used with BLD-LDAK or BLD-LDAK+Alpha models
        # NB2: thinning only required for LDAK-Thin model
        """
        ldak --calc-tagging {params.out_stem} --bfile {params.in_stem} --partition-number 4 --partition-prefix {params.partition_stem} --ignore-weights YES --window-kb 1000 --power -0.25 --max-threads {threads} > {log.log_file}
        """

use rule estimate_h2_with_human_default as estimate_h2_with_human_default_with_ighkl_annotations_for_iga_meta with:
    input:
        gwas = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/procd.assoc",
        tagging_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion}/ighkl_annotations/merged.tagging"
    output:
        multiext("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_annotations/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = "results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}/{variant_set}/{variant_type}/{ighkl_inclusion,with_ighkl}/ighkl_annotations/sumher.log"
