# Data for the manuscript 'Large-scale GWAS meta-analysis of serum antibody levels reveals distinct genetic architectures'

The `snakemake` pipeline used to generate these data can be found [here](https://github.com/twillis209/serum_ig_pipeline).

## `md5` checksums

```
1a3f4d6ac96b4ba33155aafa3e1a460a  asthma.tsv.gz
5588fde51d7f29d35fb955b81ea76e8a  celiac.tsv.gz
dad108914f27b7045f15240ef57edb50  crohns.tsv.gz
5a4ce54f92dd80a8abe74fde7fe1b823  dennis-iga.tsv.gz
be799535c1d42d2717bb356c6f326cf0  dennis-igg.tsv.gz
93b0c1bb3c5b3d04c1314713609935af  derm.tsv.gz
4c004df183ada2dd3afb5db683111c26  eldjarn-iga.tsv.gz
39fbdadbbd41d29c053f469b4bd5f5bb  eldjarn-igg.tsv.gz
8a290fd6b8fafb0ad756f8b159ec2b5b  eldjarn-igm.tsv.gz
0ec8581e5b1ea36d9dbfbd374411e098  epic-iga.tsv.gz
f2ae94fa49601845415e3a358bcc1f3b  epic-igg.tsv.gz
fddee4207aab211b7840cbed2bd7be8c  epic-igm.tsv.gz
80ef0a17979d721f267928aa885081cb  hypothy.tsv.gz
f2d464f76d1f8d6e05e9c4873ab52bde  iga-meta.tsv.gz
3ee77dd046cbf0923ae69613b400f572  igan.tsv.gz
37874f5ed1aa74a25142a69d325fd7ad  igg-meta.tsv.gz
7235b32e7879a75a6bd81450638aef84  igm-meta.tsv.gz
91309e96c5cf4af3684667ff671f8160  liu-iga.tsv.gz
23336e9af2572c42c53db4983c237668  lymphocyte-counts.tsv.gz
31306a78a15081d81122fc0ff98e5ada  ms.tsv.gz
cd85a97103a5588037ffc4d61a7b6285  pbc.tsv.gz
474f5f33c75d8c90fb12fdf07ed9a816  pietzner-iga.tsv.gz
1894e63900a3417800b2531f39b1c12c  pietzner-igg.tsv.gz
9bd364c8b3e07f5a8628f96ce9490fe4  pietzner-igm.tsv.gz
307bcced5518a0aeb007774a7a5d5bd2  psc.tsv.gz
86c3aef2e7c00e902008bcf634c23b97  ra.tsv.gz
01eb1b779a3c4751651953f9ac371545  scepanovic-iga.tsv.gz
d3957b9854fd4ef29c947faeeebf9b1a  scepanovic-igg.tsv.gz
9758780f0c0ea7e7947c0d126030e4d5  scepanovic-igm.tsv.gz
31f01dcb1c1963948ba4e6eaf6b760e2  sle.tsv.gz
7c3bb07460f19aa5c4b0b94c8955dc04  t1d.tsv.gz
a4865003b42c9924c79ea668e2e81d0e  uc.tsv.gz
```

## Column key

The GWAS summary statistics files should contain at least the following columns:

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

Note that we published the EPIC serum IgA GWAS data in an [earlier paper](https://www.sciencedirect.com/science/article/pii/S1521661624004650) but that version is superseded by this one, which has been processed with the EBI's summary statistics harmonisation pipeline rather than the one we rolled ourselves.
