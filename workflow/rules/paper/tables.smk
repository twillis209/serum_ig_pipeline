rule igh_associations_table:
    input:
        igh = "results/ig/1000kb_study_and_meta_gws_igh_lead_snps.tsv",
        igk = "results/ig/1000kb_study_and_meta_gws_igk_lead_snps.tsv"
    output:
        "results/paper/tables/igh_and_igk_associations.tsv"
    localrule: True
    script: script_path("paper/tables/igh_associations.R")

rule h2_and_rg_estimates:
    input:
        "results/ig/h2_estimates.tsv",
        "results/ig/rg_estimates.tsv"

rule iga_lead_snps:
    input:
        lead = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv",
        novel = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv",
    output:
        "results/paper/tables/iga_lead_snps.tsv"
    localrule: True
    script: script_path("paper/tables/ig_lead_snp_table.R")

use rule iga_lead_snps as igg_lead_snps with:
    input:
        lead = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv",
        novel = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv"
    output:
        "results/paper/tables/igg_lead_snps.tsv"

use rule iga_lead_snps as igm_lead_snps with:
    input:
        lead = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/1000kb_gws_annotated_lead_snps_with_study_sumstats.tsv",
        novel = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/candidate_novel_associations.tsv"
    output:
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
        "results/coloc/all_ig_pairs_with_genes.tsv"
    output:
        "results/paper/tables/ig_coloc.tsv"
    localrule: True
    run:
        daf = pd.read_csv(input[0], sep = '\t')

        daf['first_trait'] = daf['first_trait'].map(config.get('pretty_isotypes'))
        daf['second_trait'] = daf['second_trait'].map(config.get('pretty_isotypes'))

        daf.rename(columns = {'nsnps': 'No. of SNPs', "first_trait": "First isotype", "second_trait": "Second isotype", "first_snp": "First isotype's lead SNP", "second_snp": "Second isotype's lead SNP", "min_p.first": "Min. locus p-value for first isotype", "min_p.second": "Min. locus p-value for second isotype", "genes.first_snp": "Genes for first lead SNP", "genes.second_snp": "Genes for second lead SNP", "max_post": "Max. posterior hypothesis"}, inplace = True)

        daf.to_csv(output[0], sep = '\t', index = False)

rule ig_and_non_ig_coloc_results:
    input:
        "results/coloc/igg_and_asthma_results_with_genes.tsv",
        "results"
    output:
        "results/paper/tables/ig_and_non_ig_coloc.tsv"
    localrule: True
    run:
        dafs = []

        for x in input:
            daf = pd.read_csv(x, sep = '\t')

            daf['first_trait'] = daf['first_trait'].map(config.get('pretty_isotypes'))
            daf['second_trait'] = daf['second_trait'].map(config.get('gwas_datasets').get(daf['second_trait'][0]).get('pretty_phenotype'))

            daf.rename(columns = {'nsnps': 'No. of SNPs', "first_trait": "Isotype", "second_trait": "Non-Ig trait", "first_snp": "Isotype's lead SNP", "second_snp": "Non-Ig trait's lead SNP", "min_p.first": "Min. locus p-value for isotype", "min_p.second": "Min. locus p-value for non-Ig trait", "genes.first_snp": "Genes for first lead SNP", "genes.second_snp": "Genes for second lead SNP", "max_post": "Max. posterior hypothesis"}, inplace = True)

            dafs.append(daf)

        pd.concat(dafs).to_csv(output[0], sep = '\t', index = False)
