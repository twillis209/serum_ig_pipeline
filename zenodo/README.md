# Data for the manuscript 'Large-scale GWAS meta-analysis of serum antibody levels reveals distinct genetic architectures'

The `snakemake` pipeline used to generate these data can be found [here](https://github.com/twillis209/serum_ig_pipeline).

## `md5` checksums

```
0ec8581e5b1ea36d9dbfbd374411e098  epic-iga.tsv.gz
f2ae94fa49601845415e3a358bcc1f3b  epic-igg.tsv.gz
fddee4207aab211b7840cbed2bd7be8c  epic-igm.tsv.gz
f2d464f76d1f8d6e05e9c4873ab52bde  iga-meta.tsv.gz
37874f5ed1aa74a25142a69d325fd7ad  igg-meta.tsv.gz
7235b32e7879a75a6bd81450638aef84  igm-meta.tsv.gz
```

## Column key

The GWAS summary statistics files contain the following columns:

* `chromosome`: chromosome label in `hg38`
* `base_pair_location`: position in `hg38`
* `effect_allele`: effect/alternate allele
* `other_allele`: reference/other allele
* `beta`: effect estimate
* `standard_error`: standard error of effect estimate
* `p_value`: p-value
* `effect_allele_frequency`: effect allele frequency (only present in `epic-*.tsv.gz`)
* `sample_size`: maximum sample size at each SNP (only present in meta-analyses)
* `rsid`: rsID for SNP

Note that `sample_size` gives the *maximum* possible sample size at each SNP, not the actual sample size. The studies we used did not all report the number of samples/genotypes available at each SNP, so we have imputed this sample size based on the studies in which the SNP was present and the number of samples in those.

### `epic-iga.tsv.gz`

Summary statistics for the serum IgA GWAS in the EPIC-Norfolk cohort.

We published this data in an [earlier paper](https://www.sciencedirect.com/science/article/pii/S1521661624004650) but that version is superseded by this one, which has been processed with the EBI's summary statistics harmonisation pipeline rather than the one we rolled ourselves.

### `epic-igm.tsv.gz`

Summary statistics for the serum IgM GWAS in the EPIC-Norfolk cohort.

### `epic-igg.tsv.gz`

Summary statistics for the serum IgG GWAS in the EPIC-Norfolk cohort.

### `iga-meta.tsv.gz`

Summary statistics for the serum IgA GWAS meta-analysis.

### `igm-meta.tsv.gz`

Summary statistics for the serum IgM GWAS meta-analysis.

### `igg-meta.tsv.gz`

Summary statistics for the serum IgG GWAS meta-analysis.
