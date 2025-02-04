rule get_liu_lead_snps_hg19_coordinates:
    input:
        "resources/gwas/iga/liu_table_2.tsv"
    output:
        temp("results/iga_meta/existing_associations/liu_table_2_coordinates.bed")
    localrule: True
    shell: """
        tail -n +2 {input} | awk -F'\t' '{{print "chr" $2 "\t" $3 "\t" $3 "\t" $4}}' >{output}
    """

rule liftover_liu_lead_snps_to_hg38:
    input:
        liu = "results/iga_meta/existing_associations/liu_table_2_coordinates.bed",
        chain_file = "resources/liftover/hg19ToHg38.over.chain.gz"
    output:
        lifted = temp("results/iga_meta/existing_associations/liu_table_2_coordinates_lifted.bed"),
        unlifted = temp("results/iga_meta/existing_associations/liu_table_2_coordinates_unlifted.bed")
    localrule: True
    shell: "liftOver {input.liu} {input.chain_file} {output.lifted} {output.unlifted}"

rule update_liu_lead_snp_coordinates:
    input:
        bed = "results/iga_meta/existing_associations/liu_table_2_coordinates_lifted.bed",
        lead_snps = "resources/gwas/iga/liu_table_2.tsv"
    output:
        "results/iga_meta/existing_associations/liu_table_2_hg38.tsv"
    localrule: True
    run:
        bed = pd.read_csv(input.bed, sep = '\t', names = ["chrom", "start", "end", config['rsid_col']])

        lead = pd.read_csv(input.lead_snps, sep = '\t')

        merged = pd.merge(lead[['SNP', 'CHR', 'Effect Allele', 'Locus']], bed[[config['rsid_col'], 'start']], left_on = 'SNP', right_on = config['rsid_col'])

        merged = merged.rename({'CHR' : config['chr_col'], 'start': config['bp_col'], 'rsID': config['rsid_col'], 'Effect Allele': config['alt_col']}, axis = 1)

        merged[[config['chr_col'], config['bp_col'], config['rsid_col'], config['alt_col'], 'Locus']].to_csv(output[0], sep = '\t', index = False)

rule collate_existing_iga_associations:
    input:
        ebi = "resources/gwas/iga/gwas-association-downloaded_2025-02-04-EFO_0004912.tsv",
        epic = "results/harmonised_gwas/epic-iga/2000kb_gws_annotated_lead_snps.tsv",
        pietzner = "results/harmonised_gwas/pietzner-iga/2000kb_gws_annotated_lead_snps.tsv",
        gudjonsson = "results/harmonised_gwas/gudjonsson-iga/2000kb_gws_annotated_lead_snps.tsv",
        eldjarn = "results/harmonised_gwas/eldjarn-iga/2000kb_gws_annotated_lead_snps.tsv",
        scepanovic = "results/harmonised_gwas/scepanovic-iga/2000kb_gws_annotated_lead_snps.tsv",
        dennis = "results/harmonised_gwas/dennis-iga/2000kb_gws_annotated_lead_snps.tsv",
        liu = "results/iga_meta/existing_associations/liu_table_2_hg38.tsv",
        willis = "resources/gwas/iga/willis_lead_snps.tsv"
    output:
        "results/iga_meta/existing_associations.tsv"
    localrule: True
    script: script_path("iga_meta/collate_existing_associations.R")

#
#rule merge_lead_snps_with_existing_associations:
#    input:
#        existing = "results/iga_meta/existing_associations.tsv",
#        new = "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps.distance_clumped.rsIDs"
#    output:
#        "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps_and_existing_associations.tsv"
#    params:
#        window = 1e6,
#        novel = config.get('iga').get('novel')
#    localrule: True
#    script: script_path("iga_meta/merge_lead_snps_with_existing_associations.R")
#
#rule merge_all_iga_associations_with_iei_genes:
#    input:
#        iga = "results/iga_meta/with_decode/with_dennis/prescreen/gws/lead_snps_and_existing_associations.tsv",
#        iei = "resources/pid/pid_gene_coordinates.tsv"
#    output:
#        "results/iga_meta/all_associations_near_pid_genes.tsv"
#    localrule: True
#    script: script_path("iga_meta/merge_all_iga_associations_with_iei_genes.
