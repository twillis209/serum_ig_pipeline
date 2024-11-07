import re

traits_to_skip_gwas_tools = ['liu-decode-lyons-dennis-iga', 'liu-decode-lyons-iga', 'bronson-finngen-igad', 'lyons-dennis-igg', 'lim-igad']
pattern = '|'.join(re.escape(trait) for trait in traits_to_skip_gwas_tools)

wildcard_constraints:
    snp_set = 'with_mhc|sans_mhc',
    chr = "chr[0-9XY]{1,2}",
    chr_no = "[0-9XY]{1,2}",
    assembly = "hg19|hg38",
    variant_set = "all|sans_mhc|with_mhc",
    variant_type = "all|snps_only|sans_at_gc_snps_only",
    window_size = "\\d+kb",
    flank_size = "\\d+kb",
    r2 = "0_\\d+",
    download_name = '[a-z0-9\\-]+(?<!-gc)',
    input_name = fr'^(?!.*(?:{pattern}))(?!.*-gc$).*$',
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
    locus = "[\\w-]+"
