def get_rsid_and_coordinates_from_igm_lead_snps(w):
    daf = pd.read_csv(checkpoints.distance_clump_igm_meta.get(**w).output[0], sep = '\t')

    return zip(daf.rsid, daf.chromosome, daf.base_pair_location)

igm_root = "results/igm_meta/{epic_inclusion}/{scepanovic_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}"

use rule merge_iga_gwas as merge_igm_gwas with:
    input:
        epic = "results/restandardised_gwas/epic-igm.tsv.gz",
        scepanovic = "results/restandardised_gwas/scepanovic-igm.tsv.gz",
        pietzner = "resources/harmonised_gwas/pietzner-igm.tsv.gz",
        gudjonsson = "resources/harmonised_gwas/gudjonsson-igm.tsv.gz",
        eldjarn = "results/restandardised_gwas/eldjarn-igm.tsv.gz",
    output:
        "results/igm_meta/merged.tsv.gz"

use rule run_iga_meta_analysis as run_igm_meta_meta_analysis with:
    input:
        "results/igm_meta/merged.tsv.gz"
    output:
        igm_root / "meta.tsv.gz"
    params:
        isotype = 'igm'

use rule drop_selected_loci_from_iga_meta_analysis as drop_selected_loci_from_igm_meta_analysis with:
    input:
        igm_root / "meta.tsv.gz"
    output:
        igm_root / "filtered_meta.tsv.gz"

rule copy_igm_meta_to_harmonised_gwas:
    input:
        "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "resources/harmonised_gwas/igm-meta.tsv.gz"
    localrule: True
    shell: "cp {input} {output}"

checkpoint distance_clump_igm_meta:
    input:
        igm_root / "filtered_meta.tsv.gz"
    output:
        igm_root / "{window_size}_{threshold}_lead_snps.tsv"
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

use rule draw_iga_distance_clump_plot as draw_igm_distance_clump_plot with:
    input:
        igm_root / "filtered_meta.tsv.gz"
    output:
        igm_root / "{window_size}_{threshold}/{lead_rsid}_chr{chrom}_{start}_{end}.png"
    threads: 16

use rule draw_loci_from_iga_distance_clump as draw_loci_from_igm_distance_clump with:
    input:
        lambda w: [f"results/igm_meta/{{epic_inclusion}}/{{scepanovic_inclusion}}/{{pietzner_inclusion}}/{{gudjonsson_inclusion}}/{{eldjarn_inclusion}}/{{window_size}}_{{threshold}}/{rsid}_chr{chrom}_{int(int(pos)-2e6)}_{int(int(pos)+2e6)}.png" for rsid, chrom, pos in get_rsid_and_coordinates_from_igm_lead_snps(w)]
    output:
        igm_root / "{window_size}_{threshold}/locus_plots.done"

use rule collapse_clumped_iga_lead_snps as collapse_clumped_igm_lead_snps with:
    input:
        igm_root / "{window_size}_{threshold}_lead_snps.tsv"
    output:
        igm_root / "{window_size}_{threshold}_collapsed_lead_snps.tsv"
    params:
        snps_to_remove = config.get('gwas_datasets').get('igm-meta').get('lead_snps_to_remove')

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_igm_lead_snps_with_missense_and_qtl_info with:
    input:
        rules.collapse_clumped_igm_lead_snps.output
    output:
        igm_root / "{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"

use rule annotate_lead_snps_with_nearest_gene as annotate_igm_lead_snps_with_nearest_gene with:
    input:
        lead = rules.annotate_igm_lead_snps_with_missense_and_qtl_info.output,
        edb = rules.download_ensembl_db.output
    output:
        igm_root / "{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"

use rule finalise_lead_snp_annotations as finalise_igm_lead_snp_annotations with:
    input:
        rules.annotate_igm_lead_snps_with_nearest_gene
    output:
        igm_root / "{window_size}_{threshold}_annotated_lead_snps.tsv"

use rule draw_manhattan_with_lead_snp_annotation as draw_igm_meta_manhattan with:
    input:
        gwas = igm_root / "filtered_meta.tsv.gz"
    output:
        igm_root / "manhattan.png"
    params:
        title = '',
        width = 12,
        height = 4,
        ylim = [1, 1e-85]

use rule add_gnomad_queried_mafs_to_annotated_lead_snps_for_iga_meta as add_gnomad_queried_mafs_to_annotated_lead_snps_for_igm_meta with:
    input:
        igm_root / "{window_size}_{threshold}_annotated_lead_snps.tsv"
    output:
        igm_root / "{window_size}_{threshold}_annotated_lead_snps_with_gnomad_maf.tsv"

use rule add_study_sumstats_to_annotated_lead_snps_for_iga_meta as add_study_sumstats_to_annotated_lead_snps_for_igm_meta with:
    input:
        lead = igm_root / "{window_size}_{threshold}_annotated_lead_snps.tsv",
        merged = "results/igm_meta/merged.tsv.gz"
    output:
        igm_root / "{window_size}_{threshold}_annotated_lead_snps_with_study_sumstats.tsv"
    params:
        isotype = 'igm'

use rule add_novelty_flag_to_iga_lead_snps as add_novelty_flag_to_igm_lead_snps with:
    input:
        lead = igm_root / "{window_size}_{threshold}_annotated_lead_snps_with_study_sumstats.tsv",
        novel = igm_root / "candidate_novel_associations.tsv",
    output:
        igm_root / "{window_size}_{threshold}_annotated_lead_snps_with_novelty_flag.tsv"

use rule add_ieis_to_annotated_lead_snps_for_iga_meta as add_ieis_to_annotated_lead_snps_for_igm_meta with:
    input:
        lead = igm_root / "{window_size}_{threshold}_annotated_lead_snps_with_novelty_flag.tsv",
        ieis = "results/iei/ieis_by_gene.tsv"
    output:
        igm_root / "{window_size}_{threshold}_annotated_lead_snps_with_ieis.tsv"

use rule make_plink_range as make_plink_range_for_igm_meta with:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = branch(
            lambda w: w.ighkl_inclusion == 'with_ighkl',
            then = igm_root / "meta.tsv.gz",
            otherwise = igm_root / "filtered_meta.tsv.gz"
            )
    output:
        igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt"

use rule subset_reference as subset_reference_for_igm_meta with:
    input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt"
    output:
        temp(multiext(igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged", ".bed", ".bim", ".fam"))

use rule calculate_human_default_taggings as calculate_human_default_taggings_for_igm_meta with:
    input:
        multiext(igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged", ".bed", ".bim", ".fam")
    output:
        tagging_file = temp(igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged.tagging")
    log:
        log_file = igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged.tagging.log"

use rule process_sum_stats as process_sum_stats_for_igm_meta with:
    input:
        gwas_file = branch(
            lambda w: w.ighkl_inclusion == 'with_ighkl',
            then = igm_root / "meta.tsv.gz",
            otherwise = igm_root / "filtered_meta.tsv.gz"
            ),
        range_file = igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt",
    output:
        temp(igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/procd.assoc")
    params:
        N = lambda w: get_combined_sample_size_for_ig_meta(w, 'igm')

use rule estimate_h2_with_human_default as estimate_h2_with_human_default_for_igm_meta with:
    input:
        gwas = igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/procd.assoc",
        tagging_file = igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged.tagging"
    output:
        multiext(igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/sumher.", "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = igm_root / "{variant_set}/{variant_type}/{ighkl_inclusion}/sumher.log"

use rule draw_igh_locus_for_iga_datasets as draw_igh_locus_for_igm_datasets with:
    input:
        "results/igm_meta/merged.tsv.gz"
    output:
        "results/igm_meta/locuszoomr/igh.pdf"
    params:
        p_value_cols = [f"p_value.{x}" for x in config['igm_studies']],
        chrom = config['loci']['igh']['chrom'],
        start_pos = config['loci']['igh']['start'],
        stop_pos = config['loci']['igh']['stop']

use rule preprocess_sumstats_for_ldsc_munging as preprocess_igm_meta_sumstats_for_ldsc_munging with:
    input:
        sumstats = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        maf = "results/1kG/hg38/eur/snps_only/005/merged.afreq"
    output:
        temp("results/ldsc/igm/preprocessed_sumstats.tsv.gz")
