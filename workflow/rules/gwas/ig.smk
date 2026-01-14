from scipy.stats import chi2, false_discovery_control, norm

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
        expand("<iga_root>/{variant_set}/snps_only/{ighkl_inclusion}/sumher.hers", variant_set = ['sans_mhc', 'with_mhc'], ighkl_inclusion = ['with_ighkl', 'sans_ighkl']),
        expand("<igg_root>/{variant_set}/snps_only/{ighkl_inclusion}/sumher.hers", variant_set = ['sans_mhc', 'with_mhc'], ighkl_inclusion = ['with_ighkl', 'sans_ighkl']),
        expand("<igm_root>/{variant_set}/snps_only/{ighkl_inclusion}/sumher.hers", variant_set = ['sans_mhc', 'with_mhc'], ighkl_inclusion = ['with_ighkl', 'sans_ighkl'])
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
            daf['IGHKL'] = True if 'with_ighkl' in x else False
            dafs.append(daf)

            pd.concat(dafs)[['Isotype', 'Heritability model', 'MHC', 'IGHKL', 'Heritability estimate', 'Standard error']].to_csv(output[0], sep = '\t', index = False)

rule ig_imd_rg_estimates:
    input:
        cor = [f"results/ldak/ldak-thin/{x}-meta_and_{y}/inner/sans_mhc/{z}/snps_only/sumher.cors.full" for x in ["iga", "igg", "igm"] for y in config.get('imds') for z in ['with_ighkl', 'sans_ighkl']],
        log = [f"results/ldak/ldak-thin/{x}-meta_and_{y}/inner/sans_mhc/{z}/snps_only/sumher.log" for x in ["iga", "igg", "igm"] for y in config.get('imds') for z in ['with_ighkl', 'sans_ighkl']]
    output:
        "results/ig/imd_rg_estimates.tsv",
    localrule: True
    run:
        dafs = []

        for x in input.cor:
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
            daf['IGHKL'] = True if 'with_ighkl' in x else False
            dafs.append(daf)

        daf = pd.concat(dafs)

        df_wide = daf.pivot(
            index=['Isotype', 'Immune phenotype', 'Heritability model', 'MHC'],
            columns='IGHKL',
            values=['Genetic correlation estimate', 'Standard error', 'p-value']
        )

        # Flatten the MultiIndex columns (e.g., ('rg', True) -> 'rg_IGHKL_True')
        df_wide.columns = [f"{val}_IGHKL_{col}" for val, col in df_wide.columns]
        df_wide = df_wide.reset_index()

        # Calculate FDR based specifically on the IGHKL = False p-values
        # Using the flattened column name 'p-value_IGHKL_False'
        df_wide['FDR'] = false_discovery_control(df_wide['p-value_IGHKL_False'], method='bh')

        df_wide.rename(columns={
            'Genetic correlation estimate_IGHKL_True': 'Genetic correlation estimate (with IGHKL)',
            'Genetic correlation estimate_IGHKL_False': 'Genetic correlation estimate (without IGHKL)',
            'Standard error_IGHKL_True': 'Standard error (with IGHKL)',
            'Standard error_IGHKL_False': 'Standard error (without IGHKL)',
            'p-value_IGHKL_True': 'p-value (with IGHKL)',
            'p-value_IGHKL_False': 'p-value (without IGHKL)'
        }, inplace=True)

        df_wide = df_wide[['Isotype', 'Immune phenotype',
                           'Genetic correlation estimate (with IGHKL)', 'Standard error (with IGHKL)', 'p-value (with IGHKL)',
                           'Genetic correlation estimate (without IGHKL)', 'Standard error (without IGHKL)', 'p-value (without IGHKL)',
                           'FDR']]

        # Save to output
        df_wide.to_csv(output[0], sep='\t', index=False)

rule ig_rg_estimates:
    input:
        [f"results/ldak/ldak-thin/{x}/inner/{y}/{z}/snps_only/sumher.cors.full" for x in ["iga_and_igm", "iga_and_igg", "igg_and_igm"] for y in ['sans_mhc', 'with_mhc'] for z in ['with_ighkl', 'sans_ighkl']]
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
            daf['IGHKL'] = True if 'with_ighkl' in x else False
            dafs.append(daf)

        pd.concat(dafs)[['First isotype', 'Second isotype', 'Heritability model', 'MHC', 'IGHKL', 'Genetic correlation estimate', 'Standard error', 'p-value']].to_csv(output[0], sep = '\t', index = False)

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
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv"
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

rule z_test_h2_estimate_differences:
    input:
        "results/ig/h2_estimates.tsv"
    output:
        "results/ig/h2_estimate_z_tests.tsv"
    localrule: True
    run:
        df = pd.read_csv(input[0], sep = '\t')

        pivot = df.pivot_table(index=["Isotype", "MHC"], columns="IGHKL", values=["Heritability estimate", "Standard error"])

        h2_diff = pivot['Heritability estimate'][True] - pivot['Heritability estimate'][False]
        se_diff = np.sqrt(pivot['Standard error'][True]**2 + pivot['Standard error'][False]**2)
        z_score = h2_diff / se_diff
        p_values = 2 * (1 - norm.cdf(np.abs(z_score)))

        # 4. Format the final table
        results = pd.DataFrame({
            'h2_with_IGHKL': pivot['Heritability estimate'][True],
            'h2_no_IGHKL': pivot['Heritability estimate'][False],
            'Difference': h2_diff,
            'SE_diff': se_diff,
            'P_Value': p_values
        }).reset_index()

        results.to_csv(output[0], sep = '\t', index = False)
