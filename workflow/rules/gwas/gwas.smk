import pandas as pd
import os

rule join_pair_gwas:
    input:
        A = "results/processed_gwas/{trait_A}.tsv.gz",
        B = "results/processed_gwas/{trait_B}.tsv.gz"
    output:
        AB = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/merged.tsv.gz")
    threads: 8
    resources:
        runtime = 15
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
        join = lambda wildcards: wildcards.join,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P',
        beta_col = 'BETA',
        se_col = 'SE',
        id_col = 'SNPID'
    group: "gwas"
    script:
        script_path("gwas/join_pair_gwas_stats.R")

rule join_multiple_gwas_for_cfdr:
    input:
        principal_trait_gwas_file = "results/processed_gwas/{prin_trait}.tsv.gz",
        auxiliary_trait_gwas_files = lambda wildcards: [f"results/processed_gwas/{x}.tsv.gz" for x in wildcards.aux_traits.split('_and_')],
        pruned_auxiliary_trait_gwas_files = lambda wildcards: [f"results/merged_gwas/{wildcards.prin_trait}_and_{x}/inner/{wildcards.variant_set}/snps_only/{wildcards.window_size}_1_{wildcards.r2}/pruned.tsv" for x in wildcards.aux_traits.split('_and_')]
    output:
        temp("results/merged_gwas/{prin_trait}_on_{aux_traits}/left/{variant_set}/snps_only/{window_size}_1_{r2}/merged_with_prune.tsv.gz")
    threads: 8
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
        bp_col = 'BP38',
        chr_col = 'CHR38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P',
        snp_col = 'SNPID'
    resources:
        runtime = 10
    group: "gwas"
    script:
        script_path("gwas/join_multiple_gwas_stats.R")

rule make_plink_range:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = "results/processed_gwas/{trait}.tsv.gz"
    output:
        "results/processed_gwas/{trait}/{variant_set}/{variant_type}/matching_ids.txt"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
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

rule draw_qqplot:
    input:
        gwas = "results/processed_gwas/{trait}.tsv.gz",
        gif = "results/processed_gwas/gif/{trait}.tsv"
    output:
        "results/processed_gwas/plots/{trait}_qqplot.png"
    params:
        prin_col = 'P',
        pretty_trait = '', #lambda wildcards: get_metadata_field(wildcards.trait, 'pretty_name'),
        pretty_author = '' # lambda wildcards: get_metadata_field(wildcards.trait, 'First_Author')
    threads: 8
    resources:
        runtime = 10
    group: "gwas"
    script:
        script_path("gwas/plot_qqplot.R")

rule compute_genomic_inflation_factor:
    input:
        "results/processed_gwas/{trait}.tsv.gz"
    output:
        "results/processed_gwas/gif/{trait}_{variant_set,with_mhc|sans_mhc}.tsv"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = 'P',
        percentiles = [10, 20, 50, 75],
        controls = lambda w: get_metadata_field(w.trait, 'N0'),
        cases = lambda w: get_metadata_field(w.trait, 'N1'),
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    script:
        script_path("gwas/compute_genomic_inflation_factor.R")

checkpoint distance_clump_gwas:
    input:
        gwas = "results/processed_gwas/{trait}.tsv.gz",
    output:
        "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps.distance_clumped"
    params:
        mhc = lambda wildcards: False, #if wildcards.snp_set == 'sans_mhc' else True,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P',
        index_threshold = lambda wildcards: 5e-8 if wildcards.threshold == 'gws' else 1e-5,
        distance_window = 2e6,
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    script: script_path("gwas/distance_clump.R")

rule annotate_lead_snps:
    input:
        "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps.distance_clumped"
    output:
        annotations = "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps.distance_clumped.annotations",
        rsIDs = "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps.distance_clumped.rsIDs"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        snp_col = 'SNPID',
        ref_col = 'REF',
        alt_col = 'ALT'
    resources:
        runtime = 20
    localrule: True
    script:
        script_path("gwas/lead_snp_annotation.py")

rule draw_manhattan_with_lead_snp_annotation:
    input:
        gwas = "results/processed_gwas/{trait}.tsv.gz",
        rsIDs = "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps.distance_clumped.rsIDs"
    output:
        "results/processed_gwas/{trait}/{snp_set}/{threshold}/annotated_manhattan.distance_clumped.png"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38' ,
        p_col = 'P',
        snp_col = 'SNPID',
        title = '',
        width = 6,
        height = 4
    threads: 8
    resources:
        runtime = 20
    group: "gwas"
    script:
        script_path("gwas/plot_gwas_manhattan.R")

rule subset_gwas_for_lead_snp_neighbourhood:
    input:
        ld = "results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps/{snp_id}/neighbourhood.ld.gz",
        gwas = "results/processed_gwas/{trait}.tsv.gz"
    output:
        temp("results/processed_gwas/{trait}/{snp_set}/{threshold}/lead_snps/{snp_id}/neighbourhood.tsv.gz")
    params:
        lead_snp_id = lambda w: w.snp_id.replace('_', ':'),
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P',
        beta_col = 'BETA',
        se_col = 'SE',
        id_col = 'SNPID'
    threads: 8
    resources:
    group: "gwas"
    script:
        script_path("gwas/join_gwas_with_lead_snp_neighbourhood.R")

rule tabulate_merge_stats:
    input:
        A = "results/processed_gwas/{trait_A}.tsv.gz",
        B = "results/processed_gwas/{trait_B}.tsv.gz"
    output:
        AB_with_codes = "results/merged_gwas/{trait_A}_and_{trait_B}/{join,inner}/{variant_set}/merged_with_codes.tsv.gz",
        merge_stats = "results/merged_gwas/{trait_A}_and_{trait_B}/{join,inner}/{variant_set}/merge_stats.tsv"
    threads: 8
    resources:
        runtime = 15
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
        join = lambda wildcards: wildcards.join,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
    group: "gwas"
    script:
        script_path("gwas/tabulate_merge_stats.R")

rule subset_summary_statistics_about_variant:
    input:
        "results/processed_gwas/{trait}.tsv.gz"
    output:
        sum_stats = temp("results/processed_gwas/{trait}/{variant_id}/{window_size}/sum_stats.tsv.gz"),
        ids = temp("results/processed_gwas/{trait}/{variant_id}/{window_size}/ids.txt")
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        beta_cols = ['BETA'],
        se_cols = ['SE'],
        p_cols = ['P'],
        variant_id = lambda w: w.variant_id.replace('_', ':') if '_' in w.variant_id else w.variant_id,
        window = lambda w: int(w.window_size.replace('kb', '')) * 1000
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    script: script_path("ldlink/subset_sumstats_about_variant.R")

rule subset_1kGP_data_for_ld_matrix:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/processed_gwas/{trait}/{variant_id}/{window_size}/ids.txt"
    output:
        temp(multiext("results/processed_gwas/{trait}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam"))
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/processed_gwas/{trait}/{variant_id}/{window_size}/subset"
    threads: 16
    resources:
        runtime = 15
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --extract {input.ids} --make-pgen 'vzs' --out {params.out_stem}"

rule calculate_ld_for_subset_about_variant:
    input:
        multiext("results/processed_gwas/{trait}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/processed_gwas/{trait}/{variant_id}/{window_size}/subset.phased", ".vcor2", ".vcor2.vars"))
    params:
        stem = "results/processed_gwas/{trait}/{variant_id}/{window_size}/subset",
    threads: 16
    resources:
        runtime = 10
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.stem} vzs --r2-phased square --out {params.stem}"

rule merge_sumstats_with_r2:
    input:
        gwas = "results/processed_gwas/{trait}/{variant_id}/{window_size}/sum_stats.tsv.gz",
        ld = "results/processed_gwas/{trait}/{variant_id}/{window_size}/subset.phased.vcor2",
        ld_vars = "results/processed_gwas/{trait}/{variant_id}/{window_size}/subset.phased.vcor2.vars"
    output:
        temp("results/processed_gwas/{trait}/{variant_id}/{window_size}/ld_friends.tsv.gz")
    threads: 16
    resources:
        runtime = 20
    script: script_path("ldlink/merge_r2_with_sumstats.R")

rule draw_locuszoom_plot_without_r2:
    input:
        gwas = "results/processed_gwas/{trait}/{variant_id}/{window_size}/sum_stats.tsv.gz"
    output:
        "results/processed_gwas/{trait}/{variant_id}/{window_size}/locuszoom_sans_r2_{gene_track}.png"
    params:
        window = lambda w: int(w.window_size.replace('kb', '')) * 1000,
        with_genes = lambda w: True if w.gene_track == 'with_genes' else False
    container: "docker://twillis209/r-locuszoomr:latest"
    script: script_path("gwas/draw_locuszoom_plot.R")

use rule draw_locuszoom_plot_without_r2 as draw_locuszoom_plot_with_r2 with:
    input:
        gwas = "results/processed_gwas/{trait}/{variant_id}/{window_size}/ld_friends.tsv.gz"
    output:
        "results/processed_gwas/{trait}/{variant_id}/{window_size}/locuszoom_with_r2_{gene_track}.png"
