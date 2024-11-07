rule compute_ld_scores:
    input:
        multiext("results/1kG/{assembly}/{ancestry}/snps_only/merged_with_cm", ".bed", ".bim", ".fam")
    output:
        "results/ldsc/ld_scores/{assembly}/{ancestry}/merged.l2.M"
    params:
        in_stem = "results/1kG/{assembly}/{ancestry}/snps_only/merged_with_cm",
        out_stem = "results/ldsc/ld_scores/{assembly}/{ancestry}/merged"
    threads: 8
    resources:
        runtime = 150,
    container: None
    conda: env_path("ldsc.yaml")
    shell:
        """
        ldsc.py --bfile {params.in_stem} --l2 --ld-wind-cm 1 --out {params.out_stem}
        """

rule preprocess_sumstats_for_ldsc_munging:
    input:
        sumstats = "results/processed_gwas/{trait}.tsv.gz",
        maf = "results/1kG/hg38/eur/snps_only/merged.afreq"
    output:
        temp("results/ldsc/{trait}/preprocessed_sumstats.tsv.gz")
    params:
        chr_col = 'CHR38',
        bp_col = 'BP38',
        ref_col = 'REF',
        alt_col = 'ALT',
        p_col = 'P',
        beta_col = 'BETA',
        snp_col = 'SNPID'
    threads: 8
    resources:
    script:
        script_path("ldsc_and_sumher/preprocess_sumstats.R")

# TODO I think I'm misspecifying a1 and a2?
#Interpreting column names as follows:
#P:      p-Value
#ALT:    Allele 1, interpreted as ref allele for signed sumstat.
#SNPID:  Variant ID (e.g., rs number)
#REF:    Allele 2, interpreted as non-ref allele for signed sumstat.
#BETA:   Directional summary statistic as specified by --signed-sumstats.
rule munge_randomised_sum_stats:
    input:
        "results/ldsc/{trait}/preprocessed_sumstats.tsv.gz"
    output:
        temp("results/ldsc/h2/{trait}.tsv.sumstats.gz")
    params:
        output_filename = "results/ldsc/h2/{trait}.tsv",
        # NB: The '0' below gives the null value for beta
        signed_sumstats_col = "BETA,0",
        pvalue_col = lambda wildcards: "P",
        controls = lambda w: int(get_metadata_field(w.trait, 'N0')),
        cases = lambda w: int(get_metadata_field(w.trait, 'N1')),
        snp = 'SNPID',
        a1 = 'REF',
        a2 = 'ALT',
        frq = 'ALT_FREQS'
    log:
        log = "results/ldsc/h2/{trait}.tsv.log"
    threads: 1
    resources:
        runtime = 20
    container: None
    conda: env_path("ldsc.yaml")
    shell:
        """
        munge_sumstats.py --sumstats {input} --N-con {params.controls} --N-cas {params.cases} --snp {params.snp} --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 {params.a1} --a2 {params.a2} --frq {params.frq};
        """

rule estimate_h2:
    input:
        sumstats = "results/ldsc/h2/{trait}.tsv.sumstats.gz",
        ldscores = "results/ldsc/ld_scores/hg38/eur/merged.l2.M"
    output:
        "results/ldsc/h2/{trait}.log"
    params:
        ld_score_stem = "results/ldsc/ld_scores/hg38/eur/merged",
        out_stem = "results/ldsc/h2/{trait}"
    threads: 1
    resources:
        runtime = 5
    container: None
    conda: env_path("ldsc.yaml")
    shell: "ldsc.py --h2 {input.sumstats} --out {output} --ref-ld {params.ld_score_stem} --w-ld {params.ld_score_stem} --out {params.out_stem}"
