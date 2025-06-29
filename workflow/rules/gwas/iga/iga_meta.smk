def get_rsid_and_coordinates_from_iga_lead_snps(w):
    daf = pd.read_csv(checkpoints.distance_clump_iga_meta.get(**w).output[0], sep = '\t')

    return zip(daf.rsid, daf.chromosome, daf.base_pair_location)

iga_root = Path("results/iga_meta/{epic_inclusion}/{liu_inclusion}/{scepanovic_inclusion}/{dennis_inclusion}/{pietzner_inclusion}/{gudjonsson_inclusion}/{eldjarn_inclusion}")

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
    threads: 16
    resources:
        runtime = 10
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/merge_named_gwas_inputs.R")

rule run_iga_meta_analysis:
    input:
        "results/iga_meta/merged.tsv.gz"
    output:
        str(iga_root  / "meta.tsv.gz")
    params:
        isotype = 'iga'
    threads: 16
    resources:
        runtime = 30
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/run_meta_analysis.R")

rule drop_selected_loci_from_iga_meta_analysis:
    input:
        str(iga_root  / "meta.tsv.gz")
    output:
        str(iga_root  / "filtered_meta.tsv.gz")
    params:
        loci_to_drop = ['igh', 'igk', 'igl']
    threads: 16
    resources:
        runtime = 30
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/drop_loci_from_meta_analysis.R")

rule copy_iga_meta_to_harmonised_gwas:
    input:
        "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "resources/harmonised_gwas/iga-meta.tsv.gz"
    localrule: True
    shell: "cp {input} {output}"

checkpoint distance_clump_iga_meta:
    input:
        str(iga_root  / "filtered_meta.tsv.gz")
    output:
        str(iga_root  / "{window_size}_{threshold}_lead_snps.tsv")
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

rule draw_iga_distance_clump_plot:
    input:
        str(iga_root  / "filtered_meta.tsv.gz")
    output:
        str(iga_root  / "{window_size}_{threshold}/{lead_rsid}_chr{chrom}_{start}_{end}.png")
    threads: 12
    resources:
        runtime = 10
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/locuszoomr/plot_locus.R")

rule draw_loci_from_iga_distance_clump:
    input:
        lambda w: [f"results/iga_meta/{{epic_inclusion}}/{{liu_inclusion}}/{{scepanovic_inclusion}}/{{dennis_inclusion}}/{{pietzner_inclusion}}/{{gudjonsson_inclusion}}/{{eldjarn_inclusion}}/{{window_size}}_{{threshold}}/{rsid}_chr{chrom}_{int(int(pos)-2e6)}_{int(int(pos)+2e6)}.png" for rsid, chrom, pos in get_rsid_and_coordinates_from_iga_lead_snps(w)]
    output:
        str(iga_root  / "{window_size}_{threshold}/locus_plots.done")
    shell: "touch {output}"

rule collapse_clumped_iga_lead_snps:
    input:
        str(iga_root  / "{window_size}_{threshold}_lead_snps.tsv")
    output:
        str(iga_root  / "{window_size}_{threshold}_collapsed_lead_snps.tsv")
    params:
        snps_to_remove = config.get('gwas_datasets').get('iga-meta').get('lead_snps_to_remove')
    localrule: True
    run:
        pd.read_csv(input[0], sep = '\t').query("rsid not in @params.snps_to_remove").to_csv(output[0], sep = '\t', index = False)

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_iga_lead_snps_with_missense_and_qtl_info with:
    input:
        rules.collapse_clumped_iga_lead_snps.output
    output:
        str(iga_root  / "{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv")

use rule annotate_lead_snps_with_nearest_gene as annotate_iga_lead_snps_with_nearest_gene with:
    input:
        lead = rules.annotate_iga_lead_snps_with_missense_and_qtl_info.output,
        edb = "resources/gwas/ensembl_113_hsapiens_edb.sqlite"
    output:
        str(iga_root  / "{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv")

use rule finalise_lead_snp_annotations as finalise_iga_lead_snp_annotations with:
    input:
        rules.annotate_iga_lead_snps_with_nearest_gene.output
    output:
        str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps.tsv")

# NB: Taking this out for now due to timing problem with requests and redundancy of gnomad MAF estimate
rule add_gnomad_queried_mafs_to_annotated_lead_snps_for_iga_meta:
    input:
        rules.finalise_iga_lead_snp_annotations.output
    output:
        str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps_with_gnomad_maf.tsv")
    resources:
        gnomad_api_calls = 1,
        runtime = 10
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("gwas/query_lead_snp_maf_with_gnomad.py")

rule add_study_sumstats_to_annotated_lead_snps_for_iga_meta:
    input:
        lead = str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps.tsv"),
        merged = "results/iga_meta/merged.tsv.gz"
    output:
        str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps_with_study_sumstats.tsv")
    params:
        isotype = 'iga'
    threads: 8
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/add_study_sumstats_to_annotated_lead_snps.R")

rule add_novelty_flag_to_iga_lead_snps:
    input:
        lead = str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps_with_study_sumstats.tsv"),
        novel = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv",
    output:
        str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps_with_novelty_flag.tsv")
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/add_novelty_flag_to_annotated_lead_snps.R")

rule add_ieis_to_annotated_lead_snps_for_iga_meta:
    input:
        lead = str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps_with_novelty_flag.tsv"),
        ieis = "results/iei/ieis_by_gene.tsv"
    output:
        str(iga_root  / "{window_size}_{threshold}_annotated_lead_snps_with_ieis.tsv")
    params:
        flank = 2e5
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("gwas/iga_meta/add_iei_genes.R")

use rule draw_manhattan_with_lead_snp_annotation as draw_iga_meta_manhattan with:
    input:
        gwas = str(iga_root  / "filtered_meta.tsv.gz")
    output:
        str(iga_root  / "manhattan.png")
    params:
        title = '',
        width = 12,
        height = 4,
        ylim = [1, 1e-60]

use rule draw_iga_meta_manhattan as draw_annotated_iga_meta_manhattan with:
    input:
        gwas = str(iga_root  / "filtered_meta.tsv.gz"),
        lead_snps = str(iga_root  / "2000kb_gws_annotated_lead_snps.tsv")
    output:
        str(iga_root  / "annotated_manhattan.png")
    params:
        title = '',
        width = 12,
        height = 4

use rule make_plink_range as make_plink_range_for_iga_meta with:
    input:
        bim_file = "results/1kG/hg38/eur/{variant_type}/005/qc/all/merged.bim",
        gwas_file = branch(
            lambda w: w.ighkl_inclusion == 'with_ighkl',
            then = str(iga_root  / "meta.tsv.gz"),
            otherwise = str(iga_root  / "filtered_meta.tsv.gz")
        )
    output:
        str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt")

use rule subset_reference as subset_reference_for_iga_meta with:
    input:
        multiext("results/1kG/hg38/eur/{variant_type}/005/qc/all/merged", ".bed", ".bim", ".fam"),
        range_file = str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt")
    output:
        temp(multiext(str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged"), ".bed", ".bim", ".fam"))

use rule calculate_human_default_taggings as calculate_human_default_taggings_for_iga_meta with:
    input:
        multiext(str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged"), ".bed", ".bim", ".fam")
    output:
        tagging_file = temp(str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged.tagging"))
    log:
        log_file = str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged.tagging.log")

use rule process_sum_stats as process_sum_stats_for_iga_meta with:
    input:
        gwas_file = branch(
            lambda w: w.ighkl_inclusion == 'with_ighkl',
            then = str(iga_root  / "meta.tsv.gz"),
            otherwise = str(iga_root  / "filtered_meta.tsv.gz")
        ),
        range_file = str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/matching_ids.txt",)
    output:
        temp(str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/procd.assoc"))
    params:
        N = lambda w: get_combined_sample_size_for_ig_meta(w, 'iga')

use rule estimate_h2_with_human_default as estimate_h2_with_human_default_for_iga_meta with:
    input:
        gwas = str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/procd.assoc"),
        tagging_file = str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/merged.tagging")
    output:
        multiext(str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/sumher."), "cats", "cross", "enrich", "extra", "hers", "share", "taus", "progress")
    log:
        log_file = str(iga_root  / "{variant_set}/{variant_type}/{ighkl_inclusion}/sumher.log")

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

use rule preprocess_sumstats_for_ldsc_munging as preprocess_iga_meta_sumstats_for_ldsc_munging with:
    input:
        sumstats = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        maf = "results/1kG/hg38/eur/snps_only/005/merged.afreq"
    output:
        temp("results/ldsc/iga/preprocessed_sumstats.tsv.gz")
