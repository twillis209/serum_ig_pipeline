from scipy.stats import chi2
import re

# NB: LDAK-Thin model only
rule thin_predictors_for_merged_gwas:
    input:
        multiext("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output :
        thin_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/ldak/thin.in"),
        weights_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/ldak/weights.thin")
    log:
        log_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/ldak/thin.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output.thin_file, strip_suffix = '.in')
    threads: 8
    resources:
        runtime = 15
    group: "sumher"
    shell:
        """
        ldak --thin {params.out_stem} --bfile {params.in_stem} --window-prune .98 --window-kb 100 --max-threads {threads} > {log.log_file};
        awk < {output.thin_file} '{{print $1, 1}}' > {output.weights_file}
        """

rule calculate_human_default_taggings:
    input:
        multiext("results/{trait}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output :
        tagging_file = temp("results/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged.tagging"),
    log:
        log_file = "results/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged.tagging.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output.tagging_file, strip_suffix = '.tagging')
    threads: 8
    resources:
        runtime = 45
    group: "sumher"
    shell:
        # NB: weightings now only used with BLD-LDAK or BLD-LDAK+Alpha models
        # NB2: thinning only required for LDAK-Thin model
        """
        ldak --calc-tagging {params.out_stem} --bfile {params.in_stem} --ignore-weights YES --window-kb 1000 --power -0.25 --max-threads {threads} > {log.log_file}
        """

rule calculate_ldak_thin_taggings_for_merged_gwas:
    input:
        multiext("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"),
        weights_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/ldak/weights.thin"
    output:
        tagging_file = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/ldak/merged.tagging")
    log:
        log_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/ldak/merged.tagging.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.tagging')
    threads: 8
    resources:
        runtime = 15
    group: "sumher"
    shell:
        "ldak --calc-tagging {params.out_stem} --bfile {params.in_stem} --weights {input.weights_file} --window-kb 1000 --power -.25 --max-threads {threads} > {log.log_file}"

rule process_sum_stats:
    input:
        gwas_file = "results/harmonised_gwas/{trait}.tsv.gz",
        range_file = "results/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        "results/{trait}/{variant_set}/{variant_type}/merged.assoc"
    params:
        N = lambda w: config.get('gwas_datasets').get(w.trait).get('samples')
    threads: 8
    resources:
        runtime = 15,
    group: "sumher"
    script: script_path("ldsc_and_sumher/process_sum_stats.R")

rule process_sum_stats_for_merged_gwas:
    input:
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/merged.tsv.gz",
    output:
        gwas_file_A = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{trait_A}.assoc"),
        gwas_file_B = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{trait_B}.assoc")
    params:
        beta_a_col = lambda w: f'{config.get('beta_col')}.{w.trait_A}',
        beta_b_col = lambda w: f'{config.get('beta_col')}.{w.trait_B}',
        se_a_col = lambda w: f'{config.get('se_col')}.{w.trait_A}',
        se_b_col = lambda w: f'{config.get('se_col')}.{w.trait_B}',
        N_A = lambda w: config.get('gwas_datasets').get(w.trait_A).get('samples'),
        N_B = lambda w: config.get('gwas_datasets').get(w.trait_B).get('samples'),
    threads: 8
    resources:
        runtime = 15,
    group: "sumher"
    script:
        script_path("ldsc_and_sumher/process_sum_stats_for_merged_gwas.R")

rule estimate_h2_with_human_default:
    input:
        gwas = "results/{trait}/{variant_set}/{variant_type}/merged.assoc",
        tagging_file = "results/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged.tagging"
    output:
        multiext("results/ldak/human_default/{trait}/{variant_set}/{variant_type}/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = "results/ldak/human_default/{trait}/{variant_set}/{variant_type}/sumher.log"
    params:
        out_stem = subpath(output[0], strip_suffix = '.cats')
    threads: 12
    resources:
        runtime = 30
    group: "sumher"
    shell:
        "ldak --sum-hers {params.out_stem} --summary {input.gwas} --tagfile {input.tagging_file} --check-sums NO --max-threads {threads} > {log.log_file}"

rule estimate_rg_with_ldak_thin:
    input:
        tagging_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/ldak/merged.tagging",
        gwas_file_A = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{trait_A}.assoc",
        gwas_file_B = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{trait_B}.assoc"
    output:
        cors_full_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/sumher.cors.full"
    log:
        log_file = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/sumher.log"
    params:
        output_stem = subpath(output[0], strip_suffix = '.cors.full')
    threads: 8
    resources:
        runtime = 45
    group: "sumher"
    shell:
        """
        ldak --sum-cors {params.output_stem} --tagfile {input.tagging_file} --summary {input.gwas_file_A} --summary2 {input.gwas_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 --max-threads {threads} > {log.log_file}
        """

imd_trait_pairs = [f"{imd_a}_and_{imd_b}" for imd_a, imd_b in combinations(config.get('imds'), 2)]
