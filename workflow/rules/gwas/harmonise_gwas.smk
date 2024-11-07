def get_url(name):
    url = metadata_daf.loc[metadata_daf['abbrv'] == name, 'URL'].values[0]
    return url

def is_gz(name):
    return metadata_daf.loc[metadata_daf['abbrv'] == name, 'is_gz'].values[0]

def get_sha256sum(name):
    sha256 = metadata_daf.loc[metadata_daf['abbrv'] == name, 'sha256'].values[0]
    return sha256

rule download_gwas:
    output:
        ensure("resources/gwas/{download_name}.tsv.gz", sha256 = lambda wildcards: get_sha256sum(wildcards.download_name))
    params:
        url = lambda w: get_url(w.download_name),
        is_gz = lambda w: is_gz(w.download_name),
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
        output_dir = "resources/ebispot_harmoniser/reference"
    shell:
        "wget -r -P {params.output_dir} --no-parent --no-directories --continue ftp://ftp.ebi.ac.uk/pub/databases/gwas/harmonisation_resources/"
