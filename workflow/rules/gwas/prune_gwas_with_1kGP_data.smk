rule make_plink_range:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        "results/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path('gwas/prune_gwas_with_1kGP_data/make_plink_range.R')

use rule make_plink_range as make_plink_range_for_merged_gwas with:
    input:
        bim_file = "results/1kG/hg38/eur/snps_only/005/qc/all/merged.bim",
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged.tsv.gz"
    output:
        temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/matching_ids.txt")

rule subset_reference:
     input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
     output:
         temp(multiext("results/{trait}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"))
     params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.bed')
     threads: 8
     resources:
        runtime = 10
     group: "gwas"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --extract {input.range_file} --make-bed --out {params.out_stem}"

use rule subset_reference as subset_reference_for_merged_gwas with;
     input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/matching_ids.txt"
     output:
         temp(multiext("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"))
     params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.bed')

rule make_pruned_range:
    input:
        multiext("results/{trait}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output:
        temp(multiext("results/{trait}/{variant_set}/{variant_type}/{window_size}_1_{r2}/merged", ".prune.in", ".prune.out"))
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.prune.in'),
        r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
    threads: 16
    resources:
        runtime = 90
    group: "gwas"
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --indep-pairwise {wildcards.window_size} 1 {params.r2} --out {params.out_stem}"

rule prune_gwas:
     input:
        prune_file = "results/{trait}/{variant_set}/{variant_type}/{window_size}_1_{r2}/merged.prune.in",
        gwas_file = "resources/harmonised_gwas/{trait}.tsv.gz",
        maf_file = "results/1kG/hg38/eur/{variant_type}/005/merged.afreq"
     output:
         "results/{trait}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"
     threads: 16
     resources:
         runtime = 15
     group: "gwas"
     script: script_path("gwas/prune_gwas_with_1kGP_data/prune_gwas.R")
