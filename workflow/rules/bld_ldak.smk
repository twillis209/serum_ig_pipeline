rule download_bld_ldak_annotations:
    output:
        zipped = temp(ensure("resources/ldak/bld-ldak/hg19_annotations.zip", sha256 = "77cef6a817b04495c7ac811b0ad62223aa9706899707a946bbdce04bbe86ee46"))
    params:
        url = "https://genetics.ghpc.au.dk/doug/bld.zip"
    localrule: True
    shell: """wget -O {output} {params.url}"""

rule unpack_bld_ldak_annotations:
    input:
        "resources/ldak/bld-ldak/hg19_annotations.zip"
    output:
        [f"resources/ldak/bld-ldak/hg19_annotations/bld{i}" for i in range(0, 65)],
        "resources/ldak/bld-ldak/hg19_annotations/bldnames",
        "resources/ldak/bld-ldak/hg19_annotations/README.txt"
    params:
        parent = subpath(output[0], parent = True)
    localrule: True
    shell: "unzip {input} -d {params.parent}"

rule format_bld_ldak_annotations:
    input:
        "resources/ldak/bld-ldak/hg19_annotations/bld{ann_no}"
    output:
        temp("results/ldak/bld-ldak/hg19_annotations/bld{ann_no}.bed")
    localrule: True
    shell: """
        awk -F':' '{{print "chr" $1 "\t" $2 "\t" $2}}' {input} >{output}
    """

rule liftover_bld_ldak_annotations:
    input:
        bld = "results/ldak/bld-ldak/hg19_annotations/bld{ann_no}.bed",
        chain_file = "resources/liftover/hg19ToHg38.over.chain.gz"
    output:
        lifted = temp("results/ldak/bld-ldak/hg38_annotations/bld{ann_no}_lifted.bed"),
        unlifted = "results/ldak/bld-ldak/hg38_annotations/bld{ann_no}_unlifted.bed",
        renamed = "results/ldak/bld-ldak/hg38_annotations/bld{ann_no}"
    resources:
        runtime = 10
    shell: """
        liftOver {input.bld} {input.chain_file} {output.lifted} {output.unlifted};
        cp {output.lifted} {output.renamed}
    """

rule drop_pars_for_ldak_weights:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/merged", ".bim", ".bed", ".fam")
    output:
        temp(multiext("results/1kG/hg38/eur/snps_only/005/qc/all/sans_pars", ".bim", ".bed", ".fam"))
    params:
        in_stem = subpath(input[0], strip_suffix = '.bim'),
        out_stem = subpath(output[0], strip_suffix = '.bim'),
    threads: 12
    shell: "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --not-chr PAR1 PAR2 --make-bed --out {params.out_stem}"

rule cut_ldak_weights:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/sans_pars", ".bim", ".bed", ".fam")
    output:
        "results/ldak/bld-ldak/bld65/weights.predictors",
        multiext("results/ldak/bld-ldak/bld65/thin", ".dups", ".in", ".out", ".progress"),
        multiext("results/ldak/bld-ldak/bld65/section", ".details", ".number")
    params:
        in_stem = subpath(input[0], strip_suffix = '.bim'),
        out_stem = subpath(output[0], parent = True)
    threads: 12
    shell: "ldak --cut-weights {params.out_stem} --bfile {params.in_stem} --max-threads {threads}"

rule calculate_ldak_weights:
    input:
        "results/ldak/bld-ldak/bld65/weights.predictors",
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/sans_pars", ".bim", ".bed", ".fam")
    output:
        weights = "results/ldak/bld-ldak/bld65/weights.short",
        bld = "results/ldak/bld-ldak/hg38_annotations/bld65"
    params:
        stem = subpath(input[0], parent = True),
        bfile_stem = subpath(input[1], strip_suffix = ".bim")
    resources:
        runtime = 120
    threads: 12
    shell: """
        ldak --calc-weights-all {params.stem} --bfile {params.bfile_stem} --max-threads {threads};
    cp {output.weights} {output.bld}
    """

rule calculate_bld_ldak_taggings:
    input:
        multiext("results/1kG/hg38/eur/snps_only/005/qc/all/sans_pars", ".bim", ".bed", ".fam"),
        [f"results/ldak/bld-ldak/hg38_annotations/bld{x}" for x in range(1,66)]
    output:
    params:
        bfile_stem = subpath(input[0], strip_suffix = '.bim'),
        bld_stem = "results/ldak/bld-ldak/hg38_annotations/bld"
    shell: "ldak --calc-tagging BLD-LDAK --bfile {params.bfile_stem} --power -.25 --annotation-number 65 --annotation-prefix {params.bld_stem}"
