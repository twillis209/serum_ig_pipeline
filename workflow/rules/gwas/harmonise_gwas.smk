import pathlib
import yaml

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
        temp(ensure("resources/gwas/{download_name}.tsv", sha256 = lambda wildcards: config.get('gwas_datasets').get(wildcards.download_name).get('sha256')))
    params:
        url = lambda w: config.get('gwas_datasets').get(w.download_name).get('url'),
        is_gz = lambda w: config.get('gwas_datasets').get(w.download_name).get('is_gz'),
        compressed = "resources/gwas/{download_name}.tsv.gz"
    resources:
        runtime = 8
    group: "harmonise_gwas"
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
        "resources/gwas/{trait}.tsv"
    output:
        "results/gwas/gwas_ssf/{trait}.json"
    conda: env_path("gwas-sumstats-tools.yaml")
    localrule: True
    shell: "gwas-ssf format {input} --generate_config --config_out {output}"

rule write_gwas_meta_file:
    output:
        "results/gwas/gwas_ssf/{trait}.tsv-meta.yaml"
    params:
        data = lambda w: config.get('gwas_datasets').get(w.trait).get('gwas_ssf') 
    localrule: True
    run:
        with open(output[0], 'w') as fh:
            yaml.dump(params.data, fh, default_flow_style = False)

rule format_gwas:
    input:
        sumstats = "resources/gwas/{trait}.tsv",
        config = "results/gwas/gwas_ssf/{trait}.json"
    output:
        temp("results/gwas/gwas_ssf/{trait}.tsv")
    resources:
        runtime = 20
    group: "harmonise_gwas"
    conda: env_path("gwas-sumstats-tools.yaml")
    shell: "gwas-ssf format {input.sumstats} --apply_config  --config_in {input.config} -o {output}"

rule harmonise_gwas:
    input:
        config = "results/gwas/gwas_ssf/{trait}.tsv-meta.yaml",
        sumstats = "results/gwas/gwas_ssf/{trait}.tsv"
    output:
        # Bad but there are so many files to list here!
        "results/gwas/gwas_ssf/{trait}/final/{trait}.h.tsv.gz",
    params:
        launch_dir = pathlib.Path("results/gwas/gwas_ssf"),
        sumstats = lambda w: pathlib.Path(f"results/gwas/gwas_ssf/{w.trait}.tsv").resolve(),
        ref = pathlib.Path("resources/ebispot_harmoniser/reference").resolve(),
        nf_config = pathlib.Path("config/harmoniser.config").resolve(),
        version = config['ebispot_harmoniser']['version'],
        profiles = 'local,singularity'
    threads: 16
    resources:
        runtime = 90
    group: "harmonise_gwas"
    handover: True
    shell:
        """
        cd {params.launch_dir}

        nextflow \
        -c {params.nf_config} \
        run EBISPOT/gwas-sumstats-harmoniser \
        -r {params.version} \
        --ref {params.ref} \
        --harm \
        --file {params.sumstats} \
        -profile {params.profiles}
        """

rule harmonise_ig_gwas:
    input:
        expand("results/gwas/gwas_ssf/{trait}/final/{trait}.h.tsv.gz", trait = config.get('gwas_datasets'))
