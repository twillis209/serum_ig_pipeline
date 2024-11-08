rule download_human_vcf_and_rsid_reference:
    output:
        vcf = expand("resources/ebispot_harmoniser/reference/homo_sapiens-chr{chrom}.{ext}",
                     chrom = [str(x) for x in range(1,23)]+['X', 'Y', 'MT'],
                     ext = ["vcf.gz", "vcf.gz.tbi", "parquet"]
                     ),
        rsid = "resources/ebispot_harmoniser/reference/rsID.sql"
    params:
        output_dir = "resources/ebispot_harmoniser/reference",
        url = "ftp://ftp.ebi.ac.uk/pub/databases/gwas/harmonisation_resources/"
    shell:
        "wget -r -P {params.output_dir} --no-parent --no-directories --continue {params.url}"

rule download_gwas:
    output:
        ensure("resources/gwas/{download_name}.tsv", sha256 = lambda wildcards: config.get('gwas_datasets').get(wildcards.download_name).get('sha256'))
    params:
        url = lambda w: config.get('gwas_datasets').get(w.download_name).get('url'),
        is_gz = lambda w: config.get('gwas_datasets').get(w.download_name).get('is_gz'),
        uncompressed = "resources/gwas/{download_name}.tsv.gz"
    resources:
        runtime = 8
    group: "gwas"
    run:
        if params.url:
            if params.is_gz:
                shell("wget -O {params.compressed} {params.url}")
                # TODO test this
                shell("gunzip {params.compressed}")
            else:
                shell("wget -O {output} {params.url}")
        else:
            shell("exit -1")

rule generate_gwas_format_report:
    input:
        "resources/gwas/{download_name}.tsv"
    output:
        "results/gwas/gwas_ssf/{download_name}.json"
    conda: env_path("gwas-sumstats-tools.yaml")
    localrule: True
    shell: "gwas-ssf format {input} --generate_config --config_out {output}"

rule format_gwas:
    input:
        sumstats = "resources/gwas/{download_name}.tsv",
        config = "results/gwas/gwas_ssf/{download_name}.json"
    output:
        "results/gwas/gwas_ssf/{download_name}.tsv"
    resources:
        runtime = 20
    conda: env_path("gwas-sumstats-tools.yaml")
    shell: "gwas-ssf format {input.sumstats} --apply_config  --config_in {input.config} -o {output}"

