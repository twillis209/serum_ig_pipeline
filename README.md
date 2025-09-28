# `snakemake` pipeline for the manuscript 'Large-scale GWAS meta-analysis of serum antibody levels reveals distinct genetic architectures'

## Data set availability

Where possible, I've written the pipeline to download GWAS data sets from the EBI GWAS Catalog or other public URLs, but the `pietzner-*` and `eldjarn-*` data sets are not available. See the [Eldjarn](https://www.nature.com/articles/s41586-023-06563-x) and [Pietzner](https://www.science.org/doi/10.1126/science.abj1541) papers for details. In short, the Pietzner datasets can be downloaded from this [Synapse page](https://www.synapse.org/Synapse:syn51761394/wiki/622766) and you'll need to apply for access to the Eldjarn deCODE datasets then download them from AWS.

In addition, you will need to download the GWAS performed in the EPIC cohort (`epic-*`) from Zenodo. I will make these available from the GWAS Catalog eventually along with the meta-analyses, but for now both can be obtained from [Zenodo](https://zenodo.org/records/17010403).

### Harmonising raw GWAS summary statistics

The `snakemake` pipeline does not incorporate the harmonisation steps necessary to take the input GWAS from their 'raw' state to the harmonised state in which the variants have their genomic positions and alleles aligned to the `hg38` assembly (with all the attendant changes to effect estimates signs, etc). I used the [EBI's GWAS summary statistics harmoniser](https://github.com/EBISPOT/gwas-sumstats-harmoniser) to carry out these steps, but that tool uses `nextflow`, which can be integrated into `snakemake` but not in a way that I found satisfactory. I invoked the harmoniser outside the context of the `snakemake` pipeline instead. The harmonised datasets are intended to be located in `resources/harmonised_gwas`.

## How do I generate the outputs?

The tables and figures intended for publication are in the `smk` files under the `paper` directory.

## Software dependencies

I stipulate in the `Snakefile` that this pipeline requires a `snakemake` version >= 9.0. That is probably a bit strict but I like to keep up-to-date with new releases for the sake of bug fixes and new features.

All other software dependencies should be handled by `snakemake` itself, which will install them in a `conda` environment. Details of these environments can be found in the files under `workflow/envs`.
