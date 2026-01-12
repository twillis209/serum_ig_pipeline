use rule merge_iga_gwas as merge_iga_and_igg_meta with:
    input:
        iga = branch(evaluate("{ighkl_inclusion}"),
                     cases = {
                         "sans_ighkl": "<iga_root>/filtered_meta.tsv.gz",
                         "with_ighkl": "<iga_root>/meta.tsv.gz"
                         }
                     ),
        igg = branch(evaluate("{ighkl_inclusion}"),
                     cases = {
                         "sans_ighkl": "<igg_root>/filtered_meta.tsv.gz",
                         "with_ighkl": "<igg_root>/meta.tsv.gz"
                         }
                     ),
    output:
        "results/merged_gwas/iga_and_igg/{join,inner}/{variant_set}/{ighkl_inclusion}/merged.tsv.gz"
    params:
        include_sample_size = True

use rule merge_iga_gwas as merge_iga_and_igm_meta with:
    input:
        iga = branch(evaluate("{ighkl_inclusion}"),
                     cases = {
                         "sans_ighkl": "<iga_root>/filtered_meta.tsv.gz",
                         "with_ighkl": "<iga_root>/meta.tsv.gz"
                         }
                     ),
        igm = branch(evaluate("{ighkl_inclusion}"),
                     cases = {
                         "sans_ighkl": "<igm_root>/filtered_meta.tsv.gz",
                         "with_ighkl": "<igm_root>/meta.tsv.gz"
                         }
                     ),
    output:
        "results/merged_gwas/iga_and_igm/{join,inner}/{variant_set}/{ighkl_inclusion}/merged.tsv.gz"
    params:
        include_sample_size = True

use rule merge_iga_gwas as merge_igg_and_igm_meta with:
    input:
        igg = branch(evaluate("{ighkl_inclusion}"),
                     cases = {
                         "sans_ighkl": "<igg_root>/filtered_meta.tsv.gz",
                         "with_ighkl": "<igg_root>/meta.tsv.gz"
                         }
                     ),
        igm = branch(evaluate("{ighkl_inclusion}"),
                     cases = {
                         "sans_ighkl": "<igm_root>/filtered_meta.tsv.gz",
                         "with_ighkl": "<igm_root>/meta.tsv.gz"
                         }
                     ),
    output:
        "results/merged_gwas/igg_and_igm/{join,inner}/{variant_set}/{ighkl_inclusion}/merged.tsv.gz"
    params:
        include_sample_size = True

ruleorder: merge_iga_and_igg_meta > join_pair_gwas
ruleorder: merge_iga_and_igm_meta > join_pair_gwas
ruleorder: merge_igg_and_igm_meta > join_pair_gwas
