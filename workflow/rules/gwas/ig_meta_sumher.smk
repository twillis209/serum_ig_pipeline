use rule merge_iga_gwas as merge_iga_and_igg_meta with:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "results/merged_gwas/iga_and_igg/{join,inner}/{variant_set}/merged.tsv.gz"
    params:
        include_sample_size = True

use rule merge_iga_gwas as merge_iga_and_igm_meta with:
    input:
        iga = "results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "results/merged_gwas/iga_and_igm/{join,inner}/{variant_set}/merged.tsv.gz"
    params:
        include_sample_size = True

use rule merge_iga_gwas as merge_igg_and_igm_meta with:
    input:
        igg = "results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz",
        igm = "results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/filtered_meta.tsv.gz"
    output:
        "results/merged_gwas/igg_and_igm/{join,inner}/{variant_set}/merged.tsv.gz"
    params:
        include_sample_size = True

ruleorder: merge_iga_and_igg_meta > join_pair_gwas
ruleorder: merge_iga_and_igm_meta > join_pair_gwas
ruleorder: merge_igg_and_igm_meta > join_pair_gwas
