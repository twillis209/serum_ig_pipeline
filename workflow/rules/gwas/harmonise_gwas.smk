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
        "results/gwas/gwas_ssf/{trait}/{trait}.json"
    conda: env_path("gwas-sumstats-tools.yaml")
    localrule: True
    shell: "gwas-ssf format {input} --generate_config --config_out {output}"

rule write_gwas_meta_file:
    output:
        "results/gwas/gwas_ssf/{trait}/{trait}.tsv-meta.yaml"
    params:
        data = lambda w: config.get('gwas_datasets').get(w.trait).get('gwas_ssf') 
    localrule: True
    run:
        with open(output[0], 'w') as fh:
            yaml.dump(params.data, fh, default_flow_style = False)

rule format_gwas:
    input:
        sumstats = "resources/gwas/{trait}.tsv",
        config = "results/gwas/gwas_ssf/{trait}/{trait}.json"
    output:
        temp("results/gwas/gwas_ssf/{trait}/{trait}.tsv")
    resources:
        runtime = 30
    group: "harmonise_gwas"
    conda: env_path("gwas-sumstats-tools.yaml")
    shell: "gwas-ssf format {input.sumstats} --apply_config  --config_in {input.config} -o {output}"

rule harmonise_gwas:
    input:
        config = "results/gwas/gwas_ssf/{trait}/{trait}.tsv-meta.yaml",
        sumstats = "results/gwas/gwas_ssf/{trait}/{trait}.tsv"
    output:
        # Bad but there are so many files to list here!
        "results/gwas/gwas_ssf/{trait}/final/{trait}.h.tsv.gz",
    params:
        launch_dir = lambda w: pathlib.Path(f"results/gwas/gwas_ssf/{w.trait}"),
        pipeline_path = "/rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/gwas-sumstats-harmoniser",
        sumstats = lambda w: pathlib.Path(f"results/gwas/gwas_ssf/{w.trait}/{w.trait}.tsv").resolve(),
        ref = pathlib.Path("resources/ebispot_harmoniser/reference").resolve(),
        nf_config = pathlib.Path("config/harmoniser.config").resolve(),
        version = config['ebispot_harmoniser']['version'],
        profiles = 'singularity'
    threads: 16
    resources:
        runtime = 120
    group: "harmonise_gwas"
    conda: env_path("gwas_harm.yml")
    shell:
        """
        cd {params.launch_dir}

        nextflow \
        -c {params.nf_config} \
        run {params.pipeline_path} \
        --ref {params.ref} \
        --harm \
        --file {params.sumstats} \
        -profile {params.profiles}
        """

rule prep_ig_gwas:
    input:
        expand("results/gwas/gwas_ssf/{trait}/{trait}.tsv", trait = config.get('gwas_datasets'))

rule harmonise_ig_gwas:
    input:
        expand("results/gwas/gwas_ssf/{trait}/final/{trait}.h.tsv.gz", trait = config.get('gwas_datasets'))

rule estimate_sdY:
    input:
        sumstats = "results/gwas/gwas_ssf/{download_name}/final/{download_name}.h.tsv.gz",
        maf = "results/1kG/hg38/eur/snps_only/005/merged.afreq",
        prune = "results/1kG/hg38/eur/snps_only/005/qc/all/pruned/100kb_1_0_2/merged.prune.in"
    output:
        "results/gwas/gwas_ssf/sd_check/{download_name}_sdY.tsv"
    params:
        N = lambda w: config.get('gwas_datasets').get(w.download_name).get('controls'),
        no_of_snps = 20000,
        reps = 1
    threads: 8
    # TODO finish this script
    script: script_path("harmonise_gwas/estimate_sdY.R")

# TODO should have sdY standardisation here
# TODO don't we need MAF for this? Would have to get decent estimates but whence?
# TODO are there any new
rule rescale_sumstats:
    input:
        sumstats = "results/gwas/gwas_ssf/{download_name}/final/{download_name}.h.tsv.gz",
        sdY_est = "results/gwas/gwas_ssf/{download_name}/sdY.tsv"
    output:
        "results/gwas/{download_name}.tsv.gz"
    run:
        # TODO let's try polars, data.table in Python, or something
        pass
