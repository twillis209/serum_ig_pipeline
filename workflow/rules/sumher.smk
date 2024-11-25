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
        log_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/ldak/thin.log"
    params:
        in_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged",
        out_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/ldak/thin"
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
        multiext("results/processed_gwas/{trait}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output :
        tagging_file = temp("results/processed_gwas/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged.tagging"),
    log:
        log_file = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged.tagging.log"
    params:
        in_stem = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/merged",
        out_stem = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged"
    threads: 8
    resources:
        runtime = 30
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
        in_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged",
        out_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/ldak/merged"
    threads: 8
    resources:
        runtime = 15
    group: "sumher"
    shell:
        "ldak --calc-tagging {params.out_stem} --bfile {params.in_stem} --weights {input.weights_file} --window-kb 1000 --power -.25 --max-threads {threads} > {log.log_file}"

rule process_sum_stats:
    input:
        gwas_file = "results/processed_gwas/{trait}.tsv.gz",
        range_file = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
    output:
        "results/processed_gwas/{trait}/{variant_set}/{variant_type}/merged.assoc"
    params:
        N = lambda w: int(get_metadata_field(w.trait, 'N'))
    threads: 8
    resources:
        runtime = 15,
    group: "sumher"
    script:
        script_path("ldsc_and_sumher/process_sum_stats.R")

rule process_sum_stats_for_merged_gwas:
    input:
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/merged.tsv.gz",
        metadata_file = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        gwas_file_A = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{trait_A}.assoc"),
        gwas_file_B = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{trait_B}.assoc")
    params:
        trait_A = lambda wildcards: wildcards.trait_A,
        trait_B = lambda wildcards: wildcards.trait_B,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_a_col = 'BETA.A',
        beta_b_col = 'BETA.B',
        se_a_col = 'SE.A',
        se_b_col = 'SE.B'
    threads: 8
    resources:
        runtime = 15,
    group: "sumher"
    script:
        script_path("ldsc_and_sumher/process_sum_stats_for_merged_gwas.R")

rule estimate_h2_with_human_default:
    input:
        gwas = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/merged.assoc",
        tagging_file = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/ldak/human_default/merged.tagging"
    output:
        multiext("results/ldak/human_default/{trait}/{variant_set}/{variant_type}/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = "results/ldak/human_default/{trait}/{variant_set}/{variant_type}/sumher.log"
    params:
        out_stem = "results/ldak/human_default/{trait}/{variant_set}/{variant_type}/sumher",
        prevalence = lambda w: get_metadata_field(w.trait, 'prevalence'),
        case_prop = lambda w: get_metadata_field(w.trait, 'case_prop')
    threads: 8
    resources:
        runtime = 15
    group: "sumher"
    shell:
        "ldak --sum-hers {params.out_stem} --summary {input.gwas} --tagfile {input.tagging_file} --check-sums NO --prevalence {params.prevalence} --ascertainment {params.case_prop} --max-threads {threads} > {log.log_file}"

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
        output_stem = "results/ldak/ldak-thin/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/sumher"
    threads: 8
    resources:
        runtime = 20
    group: "sumher"
    shell:
        """
        ldak --sum-cors {params.output_stem} --tagfile {input.tagging_file} --summary {input.gwas_file_A} --summary2 {input.gwas_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 --max-threads {threads} > {log.log_file}
        """

imd_trait_pairs = list(chain(*[[f"{imd_traits[i]}_and_{imd_traits[j]}" for j in range(i+1,len(imd_traits))] for i in range(len(imd_traits))]))

updated_iga_imd_pairs = [x for x in imd_trait_pairs if re.search("liu-decode-lyons-dennis-iga", x)]

rule run_sumher_on_updated_iga:
    input:
        cors = [f"results/ldak/ldak-thin/{x}/inner/sans_mhc/snps_only/sumher.cors.full" for x in updated_iga_imd_pairs]

rule run_sumher_on_imds:
    input:
        cors = [f"results/ldak/ldak-thin/{x}/inner/sans_mhc/snps_only/sumher.cors.full" for x in imd_trait_pairs]
    output:
        "results/ldak/ldak-thin/combined/{variant_set}/{variant_type}/imds.tsv"
    resources:
        runtime = 30
    localrule: True
    run:
        cors_daf = compile_sumher_files(input.cors)
        #log_daf = tally_predictors_from_log_files(input.log)
        #daf = cors_daf.merge(log_daf, left_on = 'trait.B', right_on = 'trait.B')

        cors_daf.to_csv(output[0], sep = '\t', index = False)

use rule run_sumher_on_imds as run_sumher_on_trait_and_imds with:
    input:
        cors = [f"results/ldak/ldak-thin/{{trait}}_and_{x}/inner/{{variant_set}}/snps_only/sumher.cors.full" for x in imd_traits]
    output:
        "results/ldak/ldak-thin/combined/{variant_set}/snps_only/{trait}_and_imds.tsv"
    localrule: True

rule run_pid_against_updated_iga:
    input:
        "results/ldak/ldak-thin/bronson-finngen-igad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/sumher.cors.full",
        "results/ldak/ldak-thin/10kG-finngen-li-bronson-ukb-pad_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/sumher.cors.full",
        "results/ldak/ldak-thin/10kG-finngen-li-ukb-cvid_and_liu-decode-lyons-dennis-iga/inner/sans_mhc/snps_only/sumher.cors.full"

rule draw_gps_and_sumher_rg_plot_for_trait_and_imds:
    input:
        rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/{trait}_and_imds.tsv",
        gps = "results/gps/combined/{trait}/inner/sans_mhc/snps_only/1000kb_1_0_8/3000_draws/pvalues.tsv",
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        "results/ldak/ldak-thin/combined/sans_mhc/snps_only/{trait}_and_imds.png"
    localrule: True
    script: script_path("ldsc_and_sumher/draw_sumher_rg_plot.R")

rule draw_imd_heatmap:
    input:
        rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/imds.tsv",
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv"
    output:
        "results/ldak/ldak-thin/combined/sans_mhc/snps_only/imds.png"
    localrule: True
    script: script_path("ldsc_and_sumher/draw_heatmap.R")

rule collate_h2_estimates_on_liab_scale_from_combined_results:
    input:
        metadata = "resources/gwas/metadata/metadata_all_fields.tsv",
        rg = "results/ldak/ldak-thin/combined/sans_mhc/snps_only/{trait}_and_imds.tsv"
    output:
        "results/ldak/ldak-thin/combined/sans_mhc/snps_only/{trait}_and_imds_with_h2_liab.tsv"
    localrule: True
    script: script_path("ldsc_and_sumher/compute_h2_liab_estimates.R")

