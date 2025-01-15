rule ig_manhattans:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/manhattan.png" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']]
