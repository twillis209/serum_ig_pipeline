rule run_lyons_and_liu_iga_meta_analysis:
    input:
        lyons = "results/processed_gwas/iga.tsv.gz",
        liu = "results/processed_gwas/liu-iga.tsv.gz",
        liu_decode = "results/processed_gwas/liu-decode-iga.tsv.gz",
        dennis = "results/processed_gwas/dennis-iga.tsv.gz"
    output:
        liu_lyons = "results/iga_meta/without_decode/without_dennis/meta_prescreen.tsv.gz",
        liu_decode_lyons = "results/iga_meta/with_decode/without_dennis/meta_prescreen.tsv.gz",
        liu_decode_lyons_dennis = "results/iga_meta/with_decode/with_dennis/meta_prescreen.tsv.gz",
        merged = "results/iga_meta/with_decode/with_dennis/all_sum_stats.tsv.gz"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        beta_col = 'BETA',
        se_col = 'SE',
        p_col = 'P'
    threads: 8
    resources:
        runtime = 10
    group: "gwas"
    script:
        script_path("iga_meta/run_meta_analysis.R")

rule copy_iga_meta_to_processed_gwas:
    input:
        liu_lyons = "results/iga_meta/without_decode/without_dennis/meta_prescreen.tsv.gz",
        liu_decode_lyons = "results/iga_meta/with_decode/without_dennis/meta_prescreen.tsv.gz",
        liu_decode_lyons_dennis = "results/iga_meta/with_decode/with_dennis/meta_prescreen.tsv.gz"
    output:
        liu_lyons = "results/processed_gwas/liu-lyons-iga.tsv.gz",
        liu_decode_lyons = "results/processed_gwas/liu-decode-lyons-iga.tsv.gz",
        liu_decode_lyons_dennis = "results/processed_gwas/liu-decode-lyons-dennis-iga.tsv.gz"
    localrule: True
    shell: """
        cp {input.liu_lyons} {output.liu_lyons}
        cp {input.liu_decode_lyons} {output.liu_decode_lyons}
        cp {input.liu_decode_lyons_dennis} {output.liu_decode_lyons_dennis}
    """

use rule compute_genomic_inflation_factor as compute_genomic_inflation_factors_for_iga_meta with:
    input:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/meta_prescreen.tsv.gz"
    output:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{variant_set}_gif.tsv"
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        p_col = 'P',
        percentiles = [10, 20, 50, 75],
        controls = lambda w: get_metadata_field('liu-decode-lyons-iga', 'N0') if w.decode_inclusion == 'with_decode' else get_metadata_field('liu-lyons-iga', 'N0'),
        cases = lambda w: get_metadata_field('liu-decode-lyons-iga', 'N1') if w.decode_inclusion == 'with_decode' else get_metadata_field('liu-lyons-iga', 'N1')

use rule draw_qqplot as draw_iga_meta_qqplot with:
    input:
        gwas = "results/iga_meta/{decode_inclusion}/meta_prescreen.tsv.gz",
        gif = "results/iga_meta/{decode_inclusion}/gif.tsv"
    output:
        "results/iga_meta/{decode_inclusion}/meta_prescreen_qqplot.png"

checkpoint distance_clump_iga_meta:
    input:
        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/meta_{screen}.tsv.gz"
    output:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped"
    params:
        mhc = lambda wildcards: False, #if wildcards.snp_set == 'sans_mhc' else True,
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        snp_col = 'SNPID',
        p_col = 'P',
        index_threshold = lambda wildcards: 5e-8 if wildcards.threshold == 'gws' else 1e-5,
        distance_window = 2e6,
        snps_to_ignore = ['17:7581494:G:A']
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    script: script_path("iga_meta/distance_clump.R")

use rule annotate_lead_snps as annotate_iga_meta_lead_snps with:
    input:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped"
    output:
        annotations = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.annotations",
        rsIDs = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.rsIDs"

use rule draw_manhattan_with_lead_snp_annotation as draw_iga_meta_manhattan_with_lead_snp_annotation with:
    input:
        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/meta_{screen}.tsv.gz",
        rsIDs = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.rsIDs"
    output:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/annotated_manhattan.distance_clumped.png"

rule collate_existing_iga_associations:
    input:
        ebi = "resources/gwas/iga/ebi_associations.tsv",
        pietzner = "resources/gwas/iga/pietzner_associations.tsv",
        liu = "resources/gwas/iga/liu_table_2.tsv"
    output:
        "results/iga_meta/existing_associations.tsv"
    localrule: True
    script: script_path("iga_meta/collate_existing_associations.R")

use rule subset_summary_statistics_about_variant as subset_summary_statistics_about_variant_for_iga_meta with:
    input:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/meta_prescreen.tsv.gz"
    output:
        sum_stats = temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz"),
        ids = temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ids.txt")

use rule subset_1kGP_data_for_ld_matrix as subset_1kGP_data_for_iga_meta_ld_matrix with:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".pgen", ".pvar.zst", ".psam"),
        ids = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ids.txt"
    output:
        temp(multiext("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam"))
    params:
        in_stem = "results/1kG/hg38/eur/snps_only/005/qc/all/merged",
        out_stem = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset"

use rule calculate_ld_for_subset_about_variant as calculate_ld_for_subset_about_iga_meta_variant with:
    input:
        multiext("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset", ".pgen", ".pvar.zst", ".psam")
    output:
        temp(multiext("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset.phased", ".vcor2", ".vcor2.vars"))
    params:
        stem = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset",
        variant_id = lambda w: w.variant_id.replace('_', ':')

use rule fetch_rsids_for_sum_stats_about_igad_meta_variant as fetch_rsids_for_sum_stats_about_iga_meta_variant with:
    input:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz"
    output:
        temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats_with_rsids.tsv.gz")

use rule merge_sumstats_with_r2 as merge_iga_meta_sumstats_with_r2 with:
    input:
        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz",
        ld = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset.phased.vcor2",
        ld_vars = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/subset.phased.vcor2.vars"
    output:
        temp("results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ld_friends.tsv.gz")

use rule draw_locuszoom_plot_without_r2 as draw_locuszoom_plot_without_r2_for_iga_meta with:
    input:
        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/sum_stats.tsv.gz",
        rsIDs = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.rsIDs"
    output:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/locuszoom_sans_r2_{gene_track}.png"

use rule draw_locuszoom_plot_with_r2 as draw_locuszoom_plot_with_r2_for_iga_meta with:
    input:
        gwas = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/ld_friends.tsv.gz",
        rsIDs = "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/lead_snps.distance_clumped.rsIDs"
    output:
        "results/iga_meta/{decode_inclusion}/{dennis_inclusion}/{screen}/{threshold}/{variant_id}/{window_size}/locuszoom_with_r2_{gene_track}.png"

rule plot_iga_locus_tetrad:
    input:
        "results/iga_meta/with_decode/with_dennis/all_sum_stats.tsv.gz"
    output:
        "results/iga_meta/with_decode/with_dennis/prescreen/gws/locus_plots/{locus}/index_snp.png"
    params:
        index_snp_seqname = lambda w: int(config.get('iga').get('loci').get(w.locus).get('index_snp').split(':')[0]),
        index_snp_pos = lambda w: int(config.get('iga').get('loci').get(w.locus).get('index_snp').split(':')[1]),
        flank = lambda w: int(config.get('iga').get('loci').get(w.locus).get('flank').replace('kb', '')) * 1000
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    container: "docker://twillis209/r-locuszoomr:latest"
    script: script_path("iga_meta/plot_iga_locus_tetrad.R")

rule plot_all_meta_gws_loci:
    input:
        [f"results/iga_meta/with_decode/with_dennis/prescreen/gws/locus_plots/{locus}/index_snp.png" for locus in config.get('iga').get('all')]
    output:
        "results/iga_meta/with_decode/with_dennis/prescreen/gws/locus_plots/meta_gws.done"
    localrule: True
    shell: "touch {output}"

rule fetch_hg38_coordinates_for_liu_lead_snps:
    input:
        "resources/gwas/iga/liu_table_2.tsv"
    output:
    localrule: True

rule merge_liu_decode_lead_snps_with_meta_sumstats:
    input:
        liu = "resources/gwas/iga/liu_table_2.tsv",
        merged = "results/iga_meta/with_decode/with_dennis/all_sum_stats.tsv.gz"
    output:
        "results/iga_meta/with_decode/with_dennis/liu_lead_snps_with_sum_stats.tsv.gz"
    threads: 8
    localrule: True
    script: script_path("iga_meta/merge_liu_lead_snps_with_sum_stats.R")

rule merge_lead_snps_with_existing_associations:
    input:
        existing = "results/iga_meta/existing_associations.tsv",
        new = "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps.distance_clumped.rsIDs"
    output:
        "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps_and_existing_associations.tsv"
    params:
        window = 1e6,
        novel = config.get('iga').get('novel')
    localrule: True
    script: script_path("iga_meta/merge_lead_snps_with_existing_associations.R")

rule merge_all_iga_associations_with_iei_genes:
    input:
        iga = "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps_and_existing_associations.tsv",
        iei = "resources/pid/pid_gene_coordinates.tsv"
    output:
        "results/iga_meta/all_associations_near_pid_genes.tsv"
    localrule: True
    script: script_path("iga_meta/merge_all_iga_associations_with_iei_genes.R")
