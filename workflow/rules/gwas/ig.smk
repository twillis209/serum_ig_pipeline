from scipy.stats import chi2, false_discovery_control

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

rule ig_h2_estimates:
    input:
        ["results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sans_ighkl/sumher.hers",
         "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sans_ighkl/sumher.hers",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sans_ighkl/sumher.hers",
         "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sans_ighkl/sumher.hers",
         "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/with_mhc/snps_only/sans_ighkl/sumher.hers",
         "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/sans_mhc/snps_only/sans_ighkl/sumher.hers"
         ]
    output:
        "results/ig/h2_estimates.tsv"
    localrule: True
    run:
        dafs = []

        for x in input:
            isotype = config['pretty_isotypes'][x.split('/')[1].split('_')[0]]
            daf = pd.read_csv(x, sep = ' ')
            daf = daf.loc[daf['Component'] == 'Her_All', ['Heritability', 'Her_SD']]
            daf['Isotype'] = isotype
            daf = daf.rename({'Heritability': 'Heritability estimate', 'Her_SD': 'Standard error'}, axis = 1)
            daf['Heritability model'] = 'Human Default'
            daf['MHC'] = True if 'with_mhc' in x else False
            dafs.append(daf)

            pd.concat(dafs)[['Isotype', 'Heritability model', 'MHC', 'Heritability estimate', 'Standard error']].to_csv(output[0], sep = '\t', index = False)

rule ig_imd_rg_estimates:
    input:
        [f"results/ldak/ldak-thin/{x}-meta_and_{y}/inner/sans_mhc/snps_only/sumher.cors.full" for x in ["iga", "igg", "igm"] for y in config.get('imds')]
    output:
        "results/ig/imd_rg_estimates.tsv"
    localrule: True
    run:
        dafs = []

        for x in input:
            isotype = config['pretty_isotypes'].get(x.split('/')[3].split('_and_')[0])
            imd = config['gwas_datasets'].get(x.split('/')[3].split('_and_')[1]).get('pretty_phenotype')
            daf = pd.read_csv(x, sep = ' ')
            daf = daf.loc[daf['Category'] == 'ALL', ['Correlation', 'SD']]
            daf['Isotype'] = isotype
            daf['Immune phenotype'] = imd
            daf['p-value'] = chi2.sf((daf['Correlation']/daf['SD'])**2, df = 1)
            daf = daf.rename({'Correlation': 'Genetic correlation estimate', 'SD': 'Standard error'}, axis = 1)
            daf['Heritability model'] = 'LDAK-Thin'
            daf['MHC'] = True if 'with_mhc' in x else False
            dafs.append(daf)

        daf = pd.concat(dafs)

        daf['FDR'] = false_discovery_control(daf['p-value'], method = 'bh')

        daf = daf[['Isotype', 'Immune phenotype', 'Heritability model', 'MHC', 'Genetic correlation estimate', 'Standard error', 'p-value', 'FDR']]

        daf.to_csv(output[0], sep = '\t', index = False)

rule ig_rg_estimates:
    input:
        [f"results/ldak/ldak-thin/{x}/inner/{y}/snps_only/sumher.cors.full" for x in ["iga_and_igm", "iga_and_igg", "igg_and_igm"] for y in ['sans_mhc', 'with_mhc']]
    output:
        "results/ig/rg_estimates.tsv"
    localrule: True
    run:
        dafs = []

        for x in input:
            isotype_a = config['pretty_isotypes'][x.split('/')[3].split('_and_')[0]]
            isotype_b = config['pretty_isotypes'][x.split('/')[3].split('_and_')[1]]
            daf = pd.read_csv(x, sep = ' ')
            daf = daf.loc[daf['Category'] == 'ALL', ['Correlation', 'SD']]
            daf['First isotype'] = isotype_a
            daf['Second isotype'] = isotype_b
            daf['p-value'] = chi2.sf((daf['Correlation']/daf['SD'])**2, df = 1)
            daf = daf.rename({'Correlation': 'Genetic correlation estimate', 'SD': 'Standard error'}, axis = 1)
            daf['Heritability model'] = 'LDAK-Thin'
            daf['MHC'] = True if 'with_mhc' in x else False
            dafs.append(daf)

        pd.concat(dafs)[['First isotype', 'Second isotype', 'Heritability model', 'MHC', 'Genetic correlation estimate', 'Standard error', 'p-value']].to_csv(output[0], sep = '\t', index = False)

rule compile_igh_gws_hits:
    input:
         "results/harmonised_gwas/{trait}/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/harmonised_gwas/{trait}/1000kb_gws_igh_lead_snps.tsv"
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

rule compile_igh_meta_gws_hits:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps.tsv"
    output:
        "results/ig/1000kb_meta_gws_igh_lead_snps.tsv"
    params:
        chrom = config['loci']['igh']['chrom'],
        start = config['loci']['igh']['start'],
        stop = config['loci']['igh']['stop']
    localrule: True
    run:
        dafs = []

        for x in input:
            daf = pd.read_csv(x, sep = '\t')

            daf = daf.loc[
                (daf['chromosome'] == params.chrom) &
                (daf['base_pair_location'] >= params.start) &
                (daf['base_pair_location'] <= params.stop)
            ]

            daf['dataset'] = x.split('/')[1].split('_')[0]

            dafs.append(daf)

            pd.concat(dafs).to_csv(output[0], sep = '\t', index = False)

rule compile_igh_gws_hits_across_datasets:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/1000kb_gws_igh_lead_snps.tsv" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']],
        "results/ig/1000kb_meta_gws_igh_lead_snps.tsv"
    output:
        "results/ig/1000kb_study_and_meta_gws_igh_lead_snps.tsv"
    localrule: True
    run:
        pd.concat([pd.read_csv(x, sep = '\t') for x in input]).to_csv(output[0], sep = '\t', index = False)

use rule compile_igh_gws_hits as compile_igk_gws_hits with:
    output:
        "results/harmonised_gwas/{trait}/1000kb_gws_igk_lead_snps.tsv"
    params:
        chrom = config['loci']['igk']['chrom'],
        start = config['loci']['igk']['start'],
        stop = config['loci']['igk']['stop']

use rule compile_igh_meta_gws_hits as compile_igk_meta_gws_hits with:
    output:
        "results/ig/1000kb_meta_gws_igk_lead_snps.tsv"
    params:
        chrom = config['loci']['igk']['chrom'],
        start = config['loci']['igk']['start'],
        stop = config['loci']['igk']['stop']

use rule compile_igh_gws_hits_across_datasets as compile_igk_gws_hits_across_datasets with:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/1000kb_gws_igk_lead_snps.tsv" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']],
        "results/ig/1000kb_meta_gws_igk_lead_snps.tsv"
    output:
        "results/ig/1000kb_study_and_meta_gws_igk_lead_snps.tsv"

use rule compile_igh_gws_hits as compile_igl_gws_hits with:
    output:
        "results/harmonised_gwas/{trait}/1000kb_gws_igl_lead_snps.tsv"
    params:
        chrom = config['loci']['igl']['chrom'],
        start = config['loci']['igl']['start'],
        stop = config['loci']['igl']['stop']

use rule compile_igh_meta_gws_hits as compile_igl_meta_gws_hits with:
    output:
        "results/ig/1000kb_meta_gws_igl_lead_snps.tsv"
    params:
        chrom = config['loci']['igl']['chrom'],
        start = config['loci']['igl']['start'],
        stop = config['loci']['igl']['stop']

use rule compile_igh_gws_hits_across_datasets as compile_igl_gws_hits_across_datasets with:
    input:
        [[f"results/harmonised_gwas/{x}-{y}/1000kb_gws_igl_lead_snps.tsv" for x in config.get(f'{y}_studies')] for y in ['iga', 'igm', 'igg']],
        "results/ig/1000kb_meta_gws_igl_lead_snps.tsv"
    output:
        "results/ig/1000kb_study_and_meta_gws_igl_lead_snps.tsv"

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

rule estimate_phenotypic_correlations_for_epic_ig:
    input:
        "resources/epic/ig_phenotypes.csv"
    output:
        "results/ig/epic_ig_phenotypic_correlations.tsv"
    params:
        seed = 143,
        reps = 1000
    localrule: True
    conda: env_path("boot.yaml")
    script: script_path("gwas/ig/estimate_phenotypic_correlations_for_epic_ig.R")
