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
        [f"results/{trait}/{{variant_set}}/{{variant_type}}/{{window_size}}_1_{{r2}}/{seed}_sdY.tsv" for trait, seed in zip(config.get('gwas_datasets'), range(10, 10+len(config.get('gwas_datasets'))))]
    output:
        all_estimates = "results/restandardised_gwas/{variant_set}/{variant_type}/{window_size}_1_{r2}/sdY_estimates.tsv",
        sumstats = "results/restandardised_gwas/{variant_set}/{variant_type}/{window_size}_1_{r2}/sdY_estimate_sumstats.tsv"
    localrule: True
    run:
        dafs = []
        for i,x in enumerate(input):
            daf = pd.read_csv(x)
            dataset = re.search(r'results/([\w-]+)/\w+/\w+/\d+kb_1_\d+_\d+/\d+_sdY.tsv', x).groups()[0]
            daf['dataset'] = dataset
            dafs.append(daf)

        pd.concat(dafs).to_csv(output.all_estimates, sep = '\t', index = False)

        daf.groupby('dataset', as_index = False)['sdY.est'].agg(['median', 'min', 'max', 'mean']).to_csv(output.sumstats, sep = '\t', index = False)

rule restandardise_beta_and_se_using_sdY_estimate:
    input:
        "resources/harmonised_gwas/{trait}.tsv.gz"
    output:
        "results/restandardised_gwas/{trait}.tsv.gz"
    params:
        # NB: rescaling the betas and SEs if sdY estimate more than 0.05 away from 1
        sdY_estimate = lambda w: config.get('gwas_datasets').get(w.trait).get('sdY_est')
    group: "gwas"
    conda: env_path('global.yaml')
    script: script_path("gwas/phenotype_standardisation/restandardise_using_sdY.R")
