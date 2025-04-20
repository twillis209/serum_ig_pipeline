def get_rsid_and_coordinates_from_lead_snps(w):
    daf = pd.read_csv(checkpoints.distance_clump_gwas.get(**w).output[0], sep = '\t')

    return zip(daf.rsid, daf.chromosome, daf.base_pair_location)

import pandas as pd
import os

rule join_pair_gwas:
    input:
        A = "resources/harmonised_gwas/{trait_A}.tsv.gz",
        B = "resources/harmonised_gwas/{trait_B}.tsv.gz"
    output:
        AB = temp("results/merged_gwas/{trait_A}_and_{trait_B}/{join}/{variant_set}/merged.tsv.gz")
    threads: 8
    resources:
        runtime = 15
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
        join = lambda wildcards: wildcards.join,
    group: "gwas"
    script:
        script_path("gwas/join_pair_gwas_stats.R")

rule draw_qqplot:
    input:
        gwas = "resources/harmonised_gwas/{trait}.tsv.gz",
        gif = "resources/harmonised_gwas/gif/{trait}.tsv"
    output:
        "results/harmonised_gwas/plots/{trait}_qqplot.png"
    params:
        prin_col = 'p_value',
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
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        "results/{trait}/{variant_set,with_mhc|sans_mhc}_gif.tsv"
    params:
        percentiles = [10, 20, 50, 75],
        controls = lambda w: config.get('gwas_datasets').get(w.trait).get('controls'),
        cases = lambda w: config.get('gwas_datasets').get(w.trait).get('cases')
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    script:
        script_path("gwas/compute_genomic_inflation_factor.R")

checkpoint distance_clump_gwas:
    input:
        gwas = "resources/harmonised_gwas/{trait}.tsv.gz",
    output:
        "results/harmonised_gwas/{trait}/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = lambda wildcards: False, #if wildcards.snp_set == 'sans_mhc' else True,
        index_threshold = lambda wildcards: 5e-8 if wildcards.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    script: script_path("gwas/distance_clump.R")

rule draw_distance_clump_plot:
    input:
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        "results/harmonised_gwas/{trait}/{window_size}_{threshold}/{lead_rsid}_chr{chrom}_{start}_{end}.png"
    threads: 12
    resources:
        runtime = 10
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/locuszoomr/plot_locus.R")

rule draw_loci_from_distance_clump:
    input:
        lambda w: [f"results/harmonised_gwas/{{trait}}/{{window_size}}_{{threshold}}/{rsid}_chr{chrom}_{int(int(pos)-2e6)}_{int(int(pos)+2e6)}.png" for rsid, chrom, pos in get_rsid_and_coordinates_from_lead_snps(w)]
    output:
        "results/harmonised_gwas/{trait}/{window_size}_{threshold}/locus_plots.done"
    shell: "touch {output}"

rule annotate_lead_snps:
    input:
        "results/harmonised_gwas/{trait}/{window_size}_{threshold}_lead_snps.tsv"
    output:
        "results/harmonised_gwas/{trait}/{window_size}_{threshold}_annotated_lead_snps.tsv"
    params:
        no_of_top_genes = 3
    resources:
        runtime = 30
    group: "gwas"
    script:
        script_path("gwas/lead_snp_annotation.py")

rule draw_manhattan_with_lead_snp_annotation:
    input:
        gwas = "resources/harmonised_gwas/{trait}.tsv.gz",
        lead_snps = "results/harmonised_gwas/{trait}/2000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/harmonised_gwas/{trait}/manhattan.png"
    params:
        title = '',
        width = 6,
        height = 4
    threads: 16
    resources:
        runtime = 20
    group: "gwas"
    conda: env_path("global.yaml")
    script:
        script_path("gwas/plot_gwas_manhattan.R")

rule subset_gwas_for_lead_snp_neighbourhood:
    input:
        ld = "results/harmonised_gwas/{trait}/{snp_set}/{threshold}/lead_snps/{snp_id}/neighbourhood.ld.gz",
        gwas = "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        temp("resources/harmonised_gwas/{trait}/{snp_set}/{threshold}/lead_snps/{snp_id}/neighbourhood.tsv.gz")
    params:
        lead_snp_id = lambda w: w.snp_id.replace('_', ':'),
    threads: 8
    resources:
    group: "gwas"
    script:
        script_path("gwas/join_gwas_with_lead_snp_neighbourhood.R")

rule tabulate_merge_stats:
    input:
        A = "resources/harmonised_gwas/{trait_A}.tsv.gz",
        B = "resources/harmonised_gwas/{trait_B}.tsv.gz"
    output:
        AB_with_codes = "results/merged_gwas/{trait_A}_and_{trait_B}/{join,inner}/{variant_set}/merged_with_codes.tsv.gz",
        merge_stats = "results/merged_gwas/{trait_A}_and_{trait_B}/{join,inner}/{variant_set}/merge_stats.tsv"
    threads: 8
    resources:
        runtime = 15
    params:
        mhc = lambda wildcards: False if wildcards.variant_set == 'sans_mhc' else True,
        join = lambda wildcards: wildcards.join,
    group: "gwas"
    script:
        script_path("gwas/tabulate_merge_stats.R")

rule subset_summary_statistics_about_variant:
    input:
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        sum_stats = temp("results/harmonised_gwas/{trait}/{variant_id}/{window_size}/sumstats.tsv.gz"),
        ids = temp("results/harmonised_gwas/{trait}/{variant_id}/{window_size}/ids.txt")
    params:
        window = lambda w: int(w.window_size.replace('kb', '')) * 1000
    threads: 8
    resources:
        runtime = 15
    group: "gwas"
    script: script_path("gwas/subset_sumstats_about_variant.R")

rule subset_1kGP_data_for_ld_matrix:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/ids.txt"
    output:
        temp(multiext("results/harmonised_gwas/{trait}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam"))
    params:
        in_stem = subpath(input[0], strip_suffix = '.pgen'),
        out_stem = subpath(output[0], strip_suffix = '.pgen')
    threads: 16
    resources:
        runtime = 15
    conda: env_path('plink.yaml')
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.in_stem} vzs --extract {input.ids} --make-pgen 'vzs' --out {params.out_stem}"

rule calculate_ld_for_subset_about_variant:
    input:
        multiext("results/harmonised_gwas/{trait}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/harmonised_gwas/{trait}/{variant_id}/{window_size}/subset.phased", ".vcor2", ".vcor2.vars"))
    params:
        stem = subpath(input[0], strip_suffix = '.pgen')
    threads: 16
    resources:
        runtime = 10
    shell:
        "plink2 --memory {resources.mem_mb} --threads {threads} --pfile {params.stem} vzs --r2-phased square --out {params.stem}"

rule merge_sumstats_with_r2:
    input:
        gwas = "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/sumstats.tsv.gz",
        ld = "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/subset.phased.vcor2",
        ld_vars = "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/subset.phased.vcor2.vars"
    output:
        temp("results/harmonised_gwas/{trait}/{variant_id}/{window_size}/ld_friends.tsv.gz")
    threads: 16
    resources:
        runtime = 20
    script: script_path("gwas/merge_r2_with_sumstats.R")

rule draw_locuszoom_plot_without_r2:
    input:
        gwas = "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/sumstats.tsv.gz"
    output:
        "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/locuszoom_sans_r2_{gene_track}.png"
    params:
        window = lambda w: int(w.window_size.replace('kb', '')) * 1000,
        with_genes = lambda w: True if w.gene_track == 'with_genes' else False
    conda: env_path('global.yaml')
    script: script_path("gwas/locuszoomr/plot_locus.R")

use rule draw_locuszoom_plot_without_r2 as draw_locuszoom_plot_with_r2 with:
    input:
        gwas = "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/ld_friends.tsv.gz"
    output:
        "results/harmonised_gwas/{trait}/{variant_id}/{window_size}/locuszoom_with_r2_{gene_track}.png"

rule draw_locus_with_r2_from_LDlink:
    input:
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        "results/harmonised_gwas/{trait}/{locus}/locus_plot.png"
    params:
        chrom = lambda w: config.get('loci').get(w.locus).get('chrom'),
        start_pos = lambda w: config.get('loci').get(w.locus).get('start'),
        stop_pos = lambda w: config.get('loci').get(w.locus).get('stop')
    conda: env_path('global.yaml')
    script: script_path("gwas/locuszoomr/plot_locus.R")
