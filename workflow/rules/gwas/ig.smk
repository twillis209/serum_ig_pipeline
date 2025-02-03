rule ig_manhattans:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/manhattan.png" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']]

rule h2_estimates:
    input:
        ["results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sumher.hers",
         "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sumher.hers",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sumher.hers",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sumher.hers",
         "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sumher.hers",
         "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sumher.hers"
         ]

rule rg_estimates:
    input:
        "results/ldak/ldak-thin/iga-meta_and_igm-meta/inner/sans_mhc/snps_only/sumher.cors.full",
        "results/ldak/ldak-thin/iga-meta_and_igg-meta/inner/sans_mhc/snps_only/sumher.cors.full",
        "results/ldak/ldak-thin/igg-meta_and_igm-meta/inner/sans_mhc/snps_only/sumher.cors.full",

rule compile_igh_gws_hits:
    input:
        "results/harmonised_gwas/{trait}/2000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/harmonised_gwas/{trait}/2000kb_gws_igh_lead_snps.tsv"
    params:
        chrom = config['loci']['igh']['chrom'],
        start = config['loci']['igh']['start'],
        stop = config['loci']['igh']['stop']
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf = daf.loc[
            (daf['chromosome'] == params.chrom) &
            (daf['base_pair_location'] >= params.start) &
            (daf['base_pair_location'] <= params.stop)
            ]

        daf['dataset'] = wildcards.trait

        daf.to_csv(output[0], sep = '\t', index = False)

rule compile_igh_gws_hits_across_datasets:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/2000kb_gws_igh_lead_snps.tsv" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']]
    output:
        "results/harmonised_gwas/gws_igh_lead_snps.tsv"
    localrule: True
    run:
        pd.concat([pd.read_csv(x, sep = '\t') for x in input]).to_csv(output[0], sep = '\t', index = False)
