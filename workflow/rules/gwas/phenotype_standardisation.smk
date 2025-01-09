rule estimate_sdY_for_gwas:
    input:
        "results/{trait}/{variant_set}/{variant_type}/{window_size}_1_{r2}/pruned.tsv"
    output:
        "results/{trait}/{variant_set}/{variant_type}/{window_size}_1_{r2}/{seed}_sdY.tsv"
    params:
        N = lambda w: config.get('gwas_datasets').get(w.trait).get('controls'),
        no_of_sample_variants = 10000,
        reps = 10,
        seed = lambda w: int(w.seed)
    group: "gwas"
    conda: env_path('global.yaml')
    script: script_path("gwas/phenotype_standardisation/estimate_sdY_for_gwas.R")

rule estimate_sdY_for_all_datasets:
    input:
        [f"results/{trait}/sans_mhc/snps_only/1000kb_1_0_2/{seed}_sdY.tsv" for trait, seed in zip(config.get('gwas_datasets'), range(10, 10+len(config.get('gwas_datasets'))))]

rule restandardise_beta_and_se_using_sdY_estimate:
    input:
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        ""
    params:
        sdY_estimate = config.get('gwas_datasets').get(w.trait).get('sdY.est')
    group: "gwas"
    conda: env_path('global.yaml')
    script: script_path("gwas/phenotype_standardisation/restandardise_using_sdY.R")
