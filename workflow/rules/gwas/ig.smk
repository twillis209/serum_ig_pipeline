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
