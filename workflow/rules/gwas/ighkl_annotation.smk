ighkl_root = Path("results/gwas/ighkl/0")

# Filter combined IGHKL regions for each isotype

rule filter_ighkl_for_study:
    input:
        ighkl_root / "combined_ighkl_regions.tsv.gz"
    output:
        ighkl_root / "{study}/filtered_ighkl_regions.tsv.gz"
    threads: 8
    group: "gwas"
    conda: env_path("global.yaml")
    script: script_path("gwas/ighkl/filter_ighkl_for_study.R")

# IgA IGHKL lead SNP annotation

rule distance_clump_ighkl:
    input:
        ighkl_root / "{study}/filtered_ighkl_regions.tsv.gz"
    output:
        ighkl_root / "{study}/{window_size}_{threshold}_lead_snps.tsv"
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

use rule annotate_lead_snps_with_missense_and_qtl_info as annotate_ighkl_lead_snps_with_missense_and_qtl_info with:
    input:
        ighkl_root / "{study}/{window_size}_{threshold}_lead_snps.tsv"
    output:
        ighkl_root / "{study}/{window_size}_{threshold}_lead_snps_with_missense_and_qtl.tsv"

use rule annotate_lead_snps_with_nearest_gene as annotate_ighkl_lead_snps_with_nearest_gene with:
    input:
        lead = rules.annotate_ighkl_lead_snps_with_missense_and_qtl_info.output,
        edb = "resources/gwas/ensembl_113_hsapiens_edb.sqlite"
    output:
        ighkl_root / "{study}/{window_size}_{threshold}_lead_snps_with_nearest_gene.tsv"
use rule finalise_lead_snp_annotations as finalise_ighkl_lead_snp_annotations with:
    input:
        rules.annotate_ighkl_lead_snps_with_nearest_gene.output
    output:
        ighkl_root / "{study}/{window_size}_{threshold}_annotated_lead_snps.tsv"

rule annotate_ighkl_lead_snps_in_relevant_studies:
    input:
        expand(ighkl_root / "{study}/5000kb_gws_annotated_lead_snps.tsv", study = ["iga_meta", "igm_meta", "igg_meta", "eldjarn_igm", "epic_iga", "pietzner_iga", "eldjarn_iga", "epic_igg", "eldjarn_igg", "pietzner_igg"])
    output:
        "results/gwas/ighkl/0/combined_annotated_lead_snps.tsv"
    localrule: True
    run:
        import pandas as pd
        from pathlib import Path

        studies = ["iga_meta", "igm_meta", "igg_meta", "eldjarn_igm", "epic_iga", "pietzner_iga", "eldjarn_iga", "epic_igg", "eldjarn_igg", "pietzner_igg"]
        dfs = []
        for study in studies:
            path = Path(f"results/gwas/ighkl/0/{study}/5000kb_gws_annotated_lead_snps.tsv")
            if path.exists():
                df = pd.read_csv(path, sep = '\t')
                df['Study'] = config.get('gwas_datasets').get(study.replace('_', '-').get('pretty_study')
                df['Isotype'] = config.get('pretty_isotypes').get(study.split('_')[0])
                dfs.append(df)
        combined_df = pd.concat(dfs, ignore_index = True)
        combined_df.to_csv(output[0], sep = '\t', index = False)
