import re

traits_not_for_download = [x for x in config.get('gwas_datasets') if not config.get('gwas_datasets').get(x).get('url')]
pattern = '|'.join(re.escape(x) for x in traits_not_for_download)

wildcard_constraints:
    trait = '[^/]+',
    snp_set = 'with_mhc|sans_mhc',
    chr = "chr[0-9XY]{1,2}",
    chr_no = "[0-9XY]{1,2}",
    assembly = "hg19|hg38",
    variant_set = "all|sans_mhc|with_mhc|sans_pars|with_ighlk",
    ighkl_inclusion = "with_ighkl|sans_ighkl",
    variant_type = "all|snps_only|sans_at_gc_snps_only",
    window_size = "\\d+kb",
    flank_size = "\\d+kb",
    r2 = "0_\\d+",
    download_name = fr'^(?!.*(?:{pattern}))(?!.*-gc$).*$',
    maf = "\\d+",
    vmiss = "\\d+",
    post_vmiss = "\\d+",
    seed = "\\d+",
    hwe = "\\d+",
    case_hwe = "\\d+",
    control_hwe = "\\d+",
    threshold = "gws|suggestive|\\d+",
    snp_id = "\\d+_\\d+_\\w+_\\w+",
    decode_inclusion = "with_decode|without_decode",
    screen = "prescreen|postscreen",
    locus = "[\\w-]+",
    non_ig = "^!(ig(g|a|m))$"
