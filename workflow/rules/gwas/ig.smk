rule ig_manhattans:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/manhattan.png" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']],
        ["results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/manhattan.png",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/manhattan.png",
         "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/manhattan.png"
         ]

rule ig_novel_hits:
    input:
        "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv",
        "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv",
        "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv"

rule h2_estimates:
    input:
        ["results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sumher.hers",
         "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sumher.hers",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sumher.hers",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sumher.hers",
         "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sumher.hers",
         "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sumher.hers"
         ]
    output:
        "results/ig/h2_estimates.tsv"
    params:
        pretty_isotypes = {"igm_meta": "IgM", "iga_meta": "IgA", "igg_meta": "IgG"}
    localrule: True
    run:
        dafs = []

        for x in input:
            isotype = params.pretty_isotypes[x.split('/')[1]]
            daf = pd.read_csv(x, sep = ' ')
            daf = daf.loc[daf['Component'] == 'Her_All', ['Heritability', 'Her_SD']]
            daf['Isotype'] = isotype
            daf = daf.rename({'Heritability': 'Heritability estimate', 'Her_SD': 'Standard error'}, axis = 1)
            daf['Heritability model'] = 'Human Default'
            daf['MHC'] = True if 'with_mhc' in x else False
            dafs.append(daf)

        pd.concat(dafs)[['Isotype', 'Heritability model', 'MHC', 'Heritability estimate', 'Standard error']].to_csv(output[0], sep = '\t', index = False)

rule rg_estimates:
    input:
        "results/ldak/ldak-thin/iga_and_igm/inner/sans_mhc/snps_only/sumher.cors.full",
        "results/ldak/ldak-thin/iga_and_igg/inner/sans_mhc/snps_only/sumher.cors.full",
        "results/ldak/ldak-thin/igg_and_igm/inner/sans_mhc/snps_only/sumher.cors.full",

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

rule plot_beta_vs_maf_for_all_isotypes:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_gnomad_maf.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_gnomad_maf.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_gnomad_maf.tsv"
    output:
        "results/ig/beta_vs_maf.png"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("gwas/ig/plot_beta_vs_maf_at_ig_lead_snps.R")
