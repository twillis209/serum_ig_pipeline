rule download_gwas:
    output:
        ensure("resources/gwas/{download_name}.tsv.gz", sha256 = lambda wildcards: config.get('gwas_datasets').get(wildcards.download_name).get('sha256'))
    params:
        url = lambda w: config.get('gwas_datasets').get(w.download_name).get('url'),
        is_gz = lambda w: config.get('gwas_datasets').get(w.download_name).get('is_gz'),
        uncompressed = "resources/gwas/{download_name}.tsv"
    resources:
        runtime = 8
    group: "gwas"
    run:
        if params.url:
            if params.is_gz:
                shell("wget -O {output} {params.url}")
            else:
                shell("wget -O {params.uncompressed} {params.url}")
                shell("gzip {params.uncompressed}")
        else:
            shell("exit -1")

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
