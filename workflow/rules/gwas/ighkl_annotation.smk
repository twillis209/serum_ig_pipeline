ighkl_root = Path("results/gwas/ighkl/0")

# Filter combined IGHKL regions for each isotype

rule filter_ighkl_for_isotype_meta:
    input:
        ighkl_root / "combined_ighkl_regions.tsv.gz"
    output:
        ighkl_root / "{isotype}/filtered_ighkl_regions.tsv.gz"
    params:
        isotype = lambda w: w.isotype
    threads: 8
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/filter_ighkl_for_isotype.R")

# IgA IGHKL lead SNP annotation

checkpoint distance_clump_iga_ighkl:
    input:
        ighkl_root / "iga/filtered_ighkl_regions.tsv.gz"
    output:
        ighkl_root / "iga/{window_size}_{threshold}_lead_snps.tsv"
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

use rule collapse_clumped_iga_lead_snps as collapse_clumped_iga_ighkl_lead_snps with:
    input:
        ighkl_root / "iga/{window_size}_{threshold}_lead_snps.tsv"
    output:
        ighkl_root / "iga/{window_size}_{threshold}_collapsed_lead_snps.tsv"
    params:
        snps_to_remove = []

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_iga_ighkl_lead_snps_with_missense_and_qtl_info with:
    input:
        rules.collapse_clumped_iga_ighkl_lead_snps.output
    output:
        ighkl_root / "iga/{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"

use rule annotate_lead_snps_with_nearest_gene as annotate_iga_ighkl_lead_snps_with_nearest_gene with:
    input:
        lead = rules.annotate_iga_ighkl_lead_snps_with_missense_and_qtl_info.output,
        edb = "resources/gwas/ensembl_113_hsapiens_edb.sqlite"
    output:
        ighkl_root / "iga/{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"

use rule finalise_lead_snp_annotations as finalise_iga_ighkl_lead_snp_annotations with:
    input:
        rules.annotate_iga_ighkl_lead_snps_with_nearest_gene.output
    output:
        ighkl_root / "iga/{window_size}_{threshold}_annotated_lead_snps.tsv"

# IgG IGHKL lead SNP annotation

checkpoint distance_clump_igg_ighkl:
    input:
        "results/gwas/ighkl/combined_ighkl_regions.tsv.gz"
    output:
        ighkl_root / "igg/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = lambda w: False,
        index_threshold = lambda w: 5e-8 if w.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
        snps_to_ignore = [],
        p_value_col = "p_value.igg_meta"
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/distance_clump.R")

use rule collapse_clumped_iga_lead_snps as collapse_clumped_igg_ighkl_lead_snps with:
    input:
        ighkl_root / "igg/{window_size}_{threshold}_lead_snps.tsv"
    output:
        ighkl_root / "igg/{window_size}_{threshold}_collapsed_lead_snps.tsv"
    params:
        snps_to_remove = []

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_igg_ighkl_lead_snps_with_missense_and_qtl_info with:
    input:
        rules.collapse_clumped_igg_ighkl_lead_snps.output
    output:
        ighkl_root / "igg/{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"

use rule annotate_lead_snps_with_nearest_gene as annotate_igg_ighkl_lead_snps_with_nearest_gene with:
    input:
        lead = rules.annotate_igg_ighkl_lead_snps_with_missense_and_qtl_info.output,
        edb = "resources/gwas/ensembl_113_hsapiens_edb.sqlite"
    output:
        ighkl_root / "igg/{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"

use rule finalise_lead_snp_annotations as finalise_igg_ighkl_lead_snp_annotations with:
    input:
        rules.annotate_igg_ighkl_lead_snps_with_nearest_gene.output
    output:
        ighkl_root / "igg/{window_size}_{threshold}_annotated_lead_snps.tsv"

# IgM IGHKL lead SNP annotation

checkpoint distance_clump_igm_ighkl:
    input:
        "results/gwas/ighkl/combined_ighkl_regions.tsv.gz"
    output:
        ighkl_root / "igm/{window_size}_{threshold}_lead_snps.tsv"
    params:
        mhc = lambda w: False,
        index_threshold = lambda w: 5e-8 if w.threshold == 'gws' else 1e-5,
        distance_window = lambda w: int(w.window_size.split('kb')[0])*1e3,
        snps_to_ignore = [],
        p_value_col = "p_value.igm_meta"
    threads: 16
    resources:
        runtime = 5
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/distance_clump.R")

use rule collapse_clumped_iga_lead_snps as collapse_clumped_igm_ighkl_lead_snps with:
    input:
        ighkl_root / "igm/{window_size}_{threshold}_lead_snps.tsv"
    output:
        ighkl_root / "igm/{window_size}_{threshold}_collapsed_lead_snps.tsv"
    params:
        snps_to_remove = []

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_igm_ighkl_lead_snps_with_missense_and_qtl_info with:
    input:
        rules.collapse_clumped_igm_ighkl_lead_snps.output
    output:
        ighkl_root / "igm/{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"

use rule annotate_lead_snps_with_nearest_gene as annotate_igm_ighkl_lead_snps_with_nearest_gene with:
    input:
        lead = rules.annotate_igm_ighkl_lead_snps_with_missense_and_qtl_info.output,
        edb = "resources/gwas/ensembl_113_hsapiens_edb.sqlite"
    output:
        ighkl_root / "igm/{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"

use rule finalise_lead_snp_annotations as finalise_igm_ighkl_lead_snp_annotations with:
    input:
        rules.annotate_igm_ighkl_lead_snps_with_nearest_gene.output
    output:
        ighkl_root / "igm/{window_size}_{threshold}_annotated_lead_snps.tsv"
