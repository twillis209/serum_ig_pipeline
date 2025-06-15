rule igh_associations_table:
    input:
        igh = "results/ig/1000kb_study_and_meta_gws_igh_lead_snps.tsv",
        igk = "results/ig/1000kb_study_and_meta_gws_igk_lead_snps.tsv"
    output:
        "results/paper/tables/igh_and_igk_associations.tsv"
    localrule: True
    conda: env_path("global.yaml")
    script: script_path("paper/tables/igh_associations.R")

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
                print(daf)
                daf['first_trait'] = daf['first_trait'].map(config.get('pretty_isotypes'))
                daf['second_trait'] = config.get('gwas_datasets').get(daf['second_trait'][0]).get('pretty_phenotype')
                daf["genes"] = daf["genes"].apply(lambda s: ",".join(sorted(set(s.split(",")))) if pd.notna(s) else s)

                daf.rename(columns = {'nsnps': 'No. of SNPs', "first_trait": "Isotype", "second_trait": "Non-Ig trait", "ig_snp": "Isotype's lead SNP", "non_ig_snp": "Non-Ig trait's lead SNP", "chromosome": "Chromosome", "ig_snp_pos": "Ig lead SNP position", "non_ig_snp_pos": "Non-Ig lead SNP position", "min_p.first": "Min. locus p-value for isotype", "min_p.second": "Min. locus p-value for non-Ig trait", "max_post": "Max. posterior hypothesis", "pearson.cor": "Pearson correlation", "ig_snp_effect_ratio": "Effect ratio at Ig lead SNP", "non_ig_snp_effect_ratio": "Effect ratio at non-Ig lead SNP"}, inplace = True)

                dafs.append(daf)

        daf = pd.concat(dafs)

        daf['Posterior odds of H4'] = daf['PP.H4.abf']/(1. - daf['PP.H4.abf'])

        daf.to_csv(output[0], sep = '\t', index = False)

rule ig_ig_coloc_hits:
    input:
        "results/paper/tables/ig_coloc.tsv"
    output:
        "results/paper/tables/ig_coloc_hits.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf[daf['First isotype'] == 'IgA']

        rsid_daf = pd.DataFrame(config['coloc']['iga_and_igm']['hits'])
        rsid_daf['First isotype'] = 'IgA'
        rsid_daf['Second isotype'] = 'IgM'

        merged = pd.merge(daf, rsid_daf, left_on = ["First isotype", "Second isotype", "First isotype's lead SNP", "Second isotype's lead SNP"], right_on = ["First isotype", "Second isotype", 'iga', 'igm'], how = 'inner')

        merged.drop(columns = ['iga', 'igm']).to_csv(output[0], sep = '\t', index = False)

rule ig_lymphocyte_counts_hits:
    input:
        "results/paper/tables/ig_and_non_ig_coloc.tsv"
    output:
        "results/paper/tables/ig_and_non_ig_coloc_hits.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf = daf[daf['Non-Ig trait'] == 'lymphocyte count']

        iga_daf = pd.DataFrame(config['coloc']['iga_and_lymphocyte-counts']['hits'])
        iga_daf.rename(columns = {'iga': 'Isotype\'s lead SNP', 'lymphocyte-counts': 'Non-Ig trait\'s lead SNP'}, inplace = True)
        iga_daf['Isotype'] = 'IgA'
        iga_daf['Non-Ig trait'] = 'lymphocyte count'
        igg_daf = pd.DataFrame(config['coloc']['igg_and_lymphocyte-counts']['hits'])
        igg_daf.rename(columns = {'igg': 'Isotype\'s lead SNP', 'lymphocyte-counts': 'Non-Ig trait\'s lead SNP'}, inplace = True)
        igg_daf['Isotype'] = 'IgG'
        igg_daf['Non-Ig trait'] = 'lymphocyte count'
        igm_daf = pd.DataFrame(config['coloc']['igm_and_lymphocyte-counts']['hits'])
        igm_daf.rename(columns = {'igm': 'Isotype\'s lead SNP', 'lymphocyte-counts': 'Non-Ig trait\'s lead SNP'}, inplace = True)
        igm_daf['Isotype'] = 'IgM'
        igm_daf['Non-Ig trait'] = 'lymphocyte count'

        rsid_daf = pd.concat([iga_daf, igg_daf, igm_daf])

        merged = pd.merge(daf, rsid_daf, on = ['Isotype', 'Non-Ig trait', 'Isotype\'s lead SNP', 'Non-Ig trait\'s lead SNP'], how = 'inner')

        merged.to_csv(output[0], sep = '\t', index = False)

rule ig_mr:
    input:
        "results/mr/igm_and_lymphocyte-counts/3000kb_lymphocyte-counts_instruments/mr.tsv",
        "results/mr/iga_and_lymphocyte-counts/3000kb_lymphocyte-counts_instruments/mr.tsv",
        "results/mr/igg_and_lymphocyte-counts/3000kb_lymphocyte-counts_instruments/mr.tsv"
