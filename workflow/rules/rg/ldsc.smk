rule write_per_chrom_bfiles_for_ld_score_estimation:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm", ".bed", ".bim", ".fam")
    output:
        temp(multiext("results/1kG/{assembly}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm/chr{chr_no}", ".bed", ".bim", ".fam"))
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.bed')
    threads: 8
    conda: env_path("plink.yaml")
    shell: "plink2 --memory {resources.mem_mb} --threads {threads} --bfile {params.in_stem} --chr {wildcards.chr_no} --make-bed --out {params.out_stem}"

rule compute_ld_scores:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/merged_with_cm/chr{chr_no}", ".bed", ".bim", ".fam")
    output:
        multiext("results/ldsc/ld_scores/{assembly}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/chr{chr_no}", ".l2.M", ".l2.ldscore.gz", ".l2.M_5_50")
    log:
        "results/ldsc/ld_scores/{assembly}/{ancestry}/{variant_type}/{maf}/qc/{variant_set}/chr{chr_no}.log"
    params:
        in_stem = subpath(input[0], strip_suffix = '.bed'),
        out_stem = subpath(output[0], strip_suffix = '.l2.M')
    threads: 8
    resources:
        runtime = 10,
    conda: env_path("ldsc.yaml")
    shell: "ldsc.py --bfile {params.in_stem} --l2 --ld-wind-cm 1 --out {params.out_stem}"

rule compute_all_ld_scores:
    input:
       [f"results/ldsc/ld_scores/hg38/eur/snps_only/005/qc/sans_pars/chr{x}.l2.M" for x in range(1,23)]

rule preprocess_sumstats_for_ldsc_munging:
    input:
        sumstats = "resources/harmonised_gwas/{trait}.tsv.gz",
        maf = "results/1kG/hg38/eur/snps_only/005/merged.afreq"
    output:
        temp("results/ldsc/{trait}/preprocessed_sumstats.tsv.gz")
    threads: 8
    conda: env_path('global.yaml')
    script: script_path("ldsc_and_sumher/preprocess_sumstats_for_ldsc.R")

# TODO I think I'm misspecifying a1 and a2?
#Interpreting column names as follows:
#P:      p-Value
#ALT:    Allele 1, interpreted as ref allele for signed sumstat.
#SNPID:  Variant ID (e.g., rs number)
#REF:    Allele 2, interpreted as non-ref allele for signed sumstat.
#BETA:   Directional summary statistic as specified by --signed-sumstats.
rule munge_randomised_sum_stats_for_continuous_trait:
    input:
        "results/ldsc/{trait}/preprocessed_sumstats.tsv.gz"
    output:
        temp("results/ldsc/{trait}/{trait}.sumstats.gz")
    log:
        temp("results/ldsc/{trait}/{trait}.log")
    params:
        output_filename = subpath(output[0], strip_suffix = '.sumstats.gz'),
        # NB: The '0' below gives the null value for beta
        signed_sumstats_col = "beta,0",
        pvalue_col = config['p_col'],
        N = lambda w: config.get('gwas_datasets').get(w.trait).get('samples'),
        snp = 'ID',
        a1 = config['ref_col'],
        a2 = config['alt_col'],
        frq = 'ALT_FREQS'
    threads: 1
    resources:
        runtime = 10
    conda: env_path("ldsc.yaml")
    shell:
        """
        munge_sumstats.py --sumstats {input} --N {params.N} --snp {params.snp} --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 {params.a1} --a2 {params.a2} --frq {params.frq};
        """

rule estimate_h2:
    input:
        sumstats = "results/ldsc/{trait}/{trait}.sumstats.gz",
        ldscores = [f"results/ldsc/ld_scores/hg38/eur/snps_only/005/qc/sans_pars/chr{x}.l2.ldscore.gz" for x in range(1,23)]
    output:
        "results/ldsc/{trait}/h2.log"
    params:
        ld_score_stem = "results/ldsc/ld_scores/hg38/eur/snps_only/005/qc/sans_pars/chr@",
        out_stem = "results/ldsc/{trait}/h2"
    threads: 1
    resources:
        runtime = 5
    conda: env_path("ldsc.yaml")
    shell: "ldsc.py --h2 {input.sumstats} --out {output} --ref-ld-chr {params.ld_score_stem} --w-ld-chr {params.ld_score_stem} --out {params.out_stem}"

rule estimate_rg:
    input:
        sumstats_A = "results/ldsc/{trait_A}/{trait_A}.sumstats.gz",
        sumstats_B = "results/ldsc/{trait_B}/{trait_B}.sumstats.gz",
        ldscores = [f"results/ldsc/ld_scores/hg38/eur/snps_only/005/qc/sans_pars/chr{x}.l2.ldscore.gz" for x in range(1,23)]
    output:
        "results/ldsc/{trait_A}_and_{trait_B}/rg.log"
    params:
        ld_score_stem = "results/ldsc/ld_scores/hg38/eur/snps_only/005/qc/sans_pars/chr@",
        out_stem = "results/ldsc/{trait_A}_and_{trait_B}/rg"
    threads: 1
    resources:
        runtime = 5
    conda: env_path("ldsc.yaml")
    shell: "ldsc.py --rg {input.sumstats_A},{input.sumstats_B} --out {output} --ref-ld-chr {params.ld_score_stem} --w-ld-chr {params.ld_score_stem} --out {params.out_stem}"
