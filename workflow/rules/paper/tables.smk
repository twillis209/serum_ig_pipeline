rule phenotype_standardisation_table:
    input:
        rules.estimate_sdY_for_all_datasets.output.sumstats
    output:
        "results/paper/tables/phenotype_standardisation.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf.rename(columns = {'median': 'Median sd(Y) estimate', 'min': 'Min. sd(Y) estimate', 'max': 'Max. sd(Y) estimate', 'mean': 'Mean sd(Y) estimate', 'std': 'sd(Y) estimate standard error'}, inplace = True)

        daf['Antibody isotype'] = daf['dataset'].apply(lambda s: s.split('-')[1]).map(config.get('pretty_isotypes'))
        daf['Phenotype transformation'] = daf['dataset'].apply(lambda s: config.get('gwas_datasets').get(s).get('phenotype_transformation'))
        daf['Study'] = daf['dataset'].apply(lambda s: config.get('gwas_datasets').get(s).get('pretty_study'))
        daf['Restandardised?'] = daf['dataset'].apply(lambda s: config.get('gwas_datasets').get(s).get('rescale'))

        # NB: includes the Gudjonsson and Liu studies we didn't use
        daf[['Study', 'Phenotype transformation', 'Antibody isotype', 'Median sd(Y) estimate', 'Min. sd(Y) estimate', 'Max. sd(Y) estimate', 'Mean sd(Y) estimate', 'sd(Y) estimate standard deviation', 'Restandardised?']].sort_values(['Study', 'Antibody isotype']).to_csv(output[0], sep = '\t', index = False)

rule igh_associations_table:
    input:
        rules.annotate_ighkl_lead_snps_in_relevant_studies.output
    output:
        "results/paper/tables/igh_and_igk_associations.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf.rename(columns = {'rsid': 'rsID',
                                'chromosome': 'Chromosome',
                                'base_pair_location': 'Position',
                                'effect_allele': 'Effect allele',
                                'other_allele': 'Other allele',
                                'sample_size': 'Sample size',
                                'genes': 'Gene(s)',
                                'nearest_gene_name': 'Nearest gene',
                                'distance_bp': 'Distance to nearest gene',
                                'missense_gene': 'Missense gene',
                                'qtl_genes': 'QTL genes',
                                'beta': 'Beta',
                                'standard_error': 'Standard error',
                                'p_value': 'p-value',
                                },
                                inplace = True)

        daf.sort_values(['Isotype', 'Chromosome'], inplace = True)

        daf[['Study', 'Isotype', 'rsID', 'Chromosome', 'Position', 'Effect allele', 'Other allele', 'Sample size', 'Gene(s)', 'Nearest gene', 'Distance to nearest gene', 'Missense gene', 'QTL genes',  'Beta', 'Standard error', 'p-value']].to_csv(output[0], sep = '\t', index = False)


rule h2_and_rg_estimates:
    input:
        "results/ig/h2_estimates.tsv",
        "results/ig/rg_estimates.tsv"

rule iga_lead_snps:
    input:
        "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_novelty_flag.tsv"
    output:
        "results/paper/tables/iga_lead_snps.tsv"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("paper/tables/ig_lead_snp_table.R")

use rule iga_lead_snps as igg_lead_snps with:
    input:
        "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_novelty_flag.tsv"
    output:
        "results/paper/tables/igg_lead_snps.tsv"

use rule iga_lead_snps as igm_lead_snps with:
    input:
        "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_novelty_flag.tsv"
    output:
        "results/paper/tables/igm_lead_snps.tsv"

rule ig_lead_snp_tables:
    input:
        "results/paper/tables/iga_lead_snps.tsv",
        "results/paper/tables/igg_lead_snps.tsv",
        "results/paper/tables/igm_lead_snps.tsv"

rule ig_novel_lead_snp_table:
    input:
        iga = "results/paper/tables/iga_lead_snps.tsv",
        igg = "results/paper/tables/igg_lead_snps.tsv",
        igm = "results/paper/tables/igm_lead_snps.tsv"
    output:
        "results/paper/tables/ig_novel_lead_snps.tsv"
    localrule: True
    run:
        iga = pd.read_csv(input.iga, sep = '\t')
        iga['Isotype'] = 'IgA'
        igg = pd.read_csv(input.igg, sep = '\t')
        igg['Isotype'] = 'IgG'
        igm = pd.read_csv(input.igm, sep = '\t')
        igm['Isotype'] = 'IgM'

        cat = pd.concat([iga, igg, igm])

        cat[cat["Novel"] == True][['Isotype', 'rsID', 'Chromosome', 'Position', 'Effect allele', 'Other allele', 'MAF', 'Sample size', 'Gene(s)', 'Nearest gene', 'Distance to nearest gene', 'Missense gene', 'QTL genes', 'Beta', 'Standard error', 'p-value', 'Q', 'Degrees of freedom', 'I2']].to_csv(output[0], sep = '\t', index = False)

rule iei_table:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_ieis.tsv",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_ieis.tsv",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_ieis.tsv"
    output:
        "results/paper/tables/iei_table.tsv"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("paper/tables/iei_table.R")

rule imd_dataset_table:
    output:
        "results/paper/tables/imd_table.tsv"
    localrule: True
    run:
        dafs = []
        for x in config['imds']:
            dafs.append(pd.DataFrame([{k: config['gwas_datasets'].get(x).get(k, None) for k in ('pretty_phenotype', 'cases', 'controls', 'samples', 'pretty_study', 'accession_no')}]))

        daf = pd.concat(dafs)

        daf.rename(columns = {'pretty_phenotype': 'Phenotype', 'cases': 'Cases', 'controls': 'Controls', 'samples': 'Samples', 'pretty_study': 'Study', 'accession_no': 'GWAS Catalog accession'}, inplace = True)

        daf.to_csv(output[0], sep = '\t', index = False)

rule ig_coloc_results:
    input:
        "results/coloc/all_ig_pairs_with_genes_and_r2.tsv"
    output:
        "results/paper/tables/ig_coloc.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf['first_trait'] = daf['first_trait'].map(config.get('pretty_isotypes'))
        daf['second_trait'] = daf['second_trait'].map(config.get('pretty_isotypes'))
        daf["genes.first_snp"] = daf["genes.first_snp"].apply(lambda s: ",".join(sorted(set(s.split(",")))) if pd.notna(s) else s)
        daf["genes.second_snp"] = daf["genes.second_snp"].apply(lambda s: ",".join(sorted(set(s.split(",")))) if pd.notna(s) else s)

        daf.rename(columns = {'nsnps': 'No. of SNPs', "first_trait": "First isotype", "second_trait": "Second isotype", "first_snp": "First isotype's lead SNP", "second_snp": "Second isotype's lead SNP", "trimmed": "Filtered", "min_p.first": "Min. locus p-value for first isotype", "min_p.second": "Min. locus p-value for second isotype", "genes.first_snp": "Genes for first lead SNP", "genes.second_snp": "Genes for second lead SNP", "max_post": "Max. posterior hypothesis", "pearson.cor": "Pearson correlation", "first_iso_lead_snp_effect_ratio": "Effect ratio at first lead SNP", "second_iso_lead_snp_effect_ratio": "Effect ratio at second lead SNP"}, inplace = True)

        daf['Posterior odds of H4'] = daf['PP.H4.abf']/(1. - daf['PP.H4.abf'])

        # No hard cut-off but this suffices to pick out those we identified as colocalised
        daf['Colocalisation'] = daf['PP.H4.abf'] >= 0.75 & daf['Filtered'] == False

        daf.to_csv(output[0], sep = '\t', index = False)

rule ig_and_non_ig_coloc_results:
    input:
        expand("results/coloc/{isotype}_and_{non_ig}/results_with_genes_and_r2.tsv",
               isotype = ["igg", "iga", "igm"],
               non_ig = config["imds"])
    output:
        "results/paper/tables/ig_and_non_ig_coloc.tsv"
    localrule: True
    run:
        dafs = []

        for x in input:
            daf = pd.read_csv(x, sep = '\t')

            if(len(daf) > 0):
                daf['first_trait'] = daf['first_trait'].map(config.get('pretty_isotypes'))
                daf['second_trait'] = config.get('gwas_datasets').get(daf['second_trait'][0]).get('pretty_phenotype')
                daf["genes"] = daf["genes"].apply(lambda s: ",".join(sorted(set(s.split(",")))) if pd.notna(s) else s)

                daf.rename(columns = {'nsnps': 'No. of SNPs', "first_trait": "Isotype", "second_trait": "Non-Ig trait", "ig_snp": "Isotype's lead SNP", "non_ig_snp": "Non-Ig trait's lead SNP","min_p.first": "Min. locus p-value for isotype", "min_p.second": "Min. locus p-value for non-Ig trait", "max_post": "Max. posterior hypothesis", "pearson.cor": "Pearson correlation", "ig_snp_effect_ratio": "Effect ratio at Ig lead SNP", "non_ig_snp_effect_ratio": "Effect ratio at non-Ig lead SNP"}, inplace = True)

                dafs.append(daf)

        daf = pd.concat(dafs)

        daf['Posterior odds of H4'] = daf['PP.H4.abf']/(1. - daf['PP.H4.abf'])

        # int_cols = ["Chromosome", "Ig lead SNP position", "Non-Ig lead SNP position"]

        # daf[int_cols] = daf[int_cols].astype("Int64")

        daf.to_csv(output[0], sep = '\t', index = False)
