rule make_plink_range:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = "results/processed_gwas/{trait}.tsv.gz"
    output:
        "results/processed_gwas/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'other_allele',
        alt_col = 'effect_allele',
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    script: script_path('gwas/make_plink_range.R')

use rule make_plink_range as make_plink_range_for_merged_gwas with:
    input:
        bim_file = "results/1kG/hg38/eur/snps_only/005/qc/all/merged.bim",
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/merged.tsv.gz"
    output:
        temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/matching_ids.txt")

rule subset_reference:
     input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
     output:
         temp(multiext("results/processed_gwas/{trait}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam"))
     params:
        in_stem = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged",
        out_stem = "results/processed_gwas/{trait}/{variant_set}/{variant_type}/merged"
     threads: 8
     resources:
        runtime = 10
     group: "gwas"
     shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --extract {input.range_file} --make-bed --out {params.out_stem}"

use rule subset_reference as subset_reference_for_merged_gwas with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/matching_ids.txt"
    output:
        temp(multiext("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/merged", ".bed", ".bim", ".fam"))
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type,snps_only}/merged"

rule make_pruned_range:
    input:
        multiext("results/processed_gwas/{trait}/merged", ".bed", ".bim", ".fam")
    output:
        temp(multiext("results/processed_gwas/{trait}/{window_size}_1_{r2}/merged", ".prune.in", ".prune.out"))
    params:
            in_stem = "results/processed_gwas/{trait}/merged",
            out_stem = "results/processed_gwas/{trait}/{window_size}_1_{r2}/merged",
            r2 = lambda wildcards: wildcards.r2.replace('_', '.'),
    threads: 16
    resources:
        runtime = 90
    group: "gwas"
    shell:
        "plink --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --indep-pairwise {wildcards.window_size} 1 {params.r2} --out {params.out_stem}"

use rule make_pruned_range as make_pruned_range_for_merged_gwas with:
    input:
        multiext("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged", ".bed", ".bim", ".fam")
    output:
        temp(multiext("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/merged", ".prune.in", ".prune.out"))
    params:
        in_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/merged",
        out_stem = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/merged",
        r2 = lambda wildcards: wildcards.r2.replace('_', '.')

rule prune_gwas:
     input:
        prune_file = "results/processed_gwas/{trait}/{window_size}_1_{r2}/merged.prune.in",
        gwas_file = "results/processed_gwas/{trait}.tsv.gz"
     output:
         "results/processed_gwas/{trait}/{window_size}_1_{r2}/pruned.tsv"
     params:
         bp_col = 'BP38',
         chr_col = 'CHR38',
         ref_col = 'REF',
         alt_col = 'ALT',
     threads: 16
     resources:
         runtime = 15
     group: "gwas"
     script: script_path("gwas/prune_gwas.R")

use rule prune_gwas as prune_merged_gwas with:
    input:
        prune_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/merged.prune.in",
        gwas_file = "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/merged.tsv.gz"
    output:
        "results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"

rule draw_manhattan:
    input:
        gwas = "results/processed_gwas/{trait}.tsv.gz"
    output:
        "results/processed_gwas/plots/{trait}_manhattan.{ext}"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = 'P',
        snp_col = 'SNPID',
        title = '',
        #title = lambda wildcards: 'GWAS of %s (%s et al.)' % (get_metadata_field(wildcards.trait, 'pretty_name'), get_metadata_field(wildcards.trait, 'First_Author')),
        width = 6,
        height = 4
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    script:
        script_path("gwas/plot_gwas_manhattan.R")
