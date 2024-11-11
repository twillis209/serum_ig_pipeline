import pathlib

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

rule write_gwas_meta_file:
    output:
    run:
        pass

rule format_gwas:
    input:
        sumstats = "resources/gwas/{download_name}.tsv",
        config = "results/gwas/gwas_ssf/{download_name}.json"
    output:
        temp("results/gwas/gwas_ssf/{download_name}.tsv")
    resources:
        runtime = 20
    conda: env_path("gwas-sumstats-tools.yaml")
    shell: "gwas-ssf format {input.sumstats} --apply_config  --config_in {input.config} -o {output}"

rule harmonise_gwas:
    input:
        config = "results/gwas/gwas_ssf/{download_name}.tsv-meta.yaml",
        sumstats = "results/gwas/gwas_ssf/{download_name}.tsv"
    output:
        # Bad but there are so many files to list here!
	"results/gwas/gwas_sff/{download_name}/final/{download_name}.h.tsv.gz",
	temp(directory("results/gwas/gwas_sff/work"))
    params:
        launch_dir = pathlib.Path("results/gwas/gwas_ssf"),
        sumstats = lambda w: pathlib.Path(f"results/gwas/gwas_ssf/{w.download_name}.tsv").resolve(),
        ref = pathlib.Path("resources/ebispot_harmoniser/reference").resolve(),
        nf_config = pathlib.Path("config/harmoniser.config").resolve(),
        version = config['ebispot_harmoniser']['version'],
        profile = 'local,singularity'
    threads: 12
    resources:
        runtime = 90
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
        -profile {params.profile}
        """

rule test_harmonise_gwas:
    output:
        "random_name/final/random_name.h.tsv.gz"
    params:
        ref = pathlib.Path("resources/ebispot_harmoniser/reference").resolve(),
        version = config['ebispot_harmoniser']['version'],
        output_dir = pathlib.Path(".").resolve()
    threads: 16
    resources:
        runtime = 120
    handover: True
    shell:
        """
        nextflow run EBISPOT/gwas-sumstats-harmoniser \
        -r {params.version} \
        --ref {params.ref} \
        --output-dir {params.output_dir} \
        -profile test,singularity
        """

