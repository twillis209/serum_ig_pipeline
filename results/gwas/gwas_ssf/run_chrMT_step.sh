coordinate_system=$(grep coordinate_system pietzner-igg.tsv-meta.yaml | awk -F ":" '{print $2}' | tr -d "[:blank:]" )
if test -z "$coordinate_system"; then coordinate="1-based"; else coordinate=$coordinate_system; fi

conda activate gwas_harm

header_args=$(/rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/gwas-sumstats-harmoniser/bin/utils.py -f /rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/gwas/gwas_ssf/pietzner-igg/1_map_to_build/MT.merged -harm_args);

harm_dir="/rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/gwas/gwas_ssf/pietzner-igg/4_harmonization"

/rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/gwas-sumstats-harmoniser/bin/main_pysam.py \
--sumstats /rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/gwas/gwas_ssf/pietzner-igg/1_map_to_build/MT.merged \
--vcf /rds/project/rds-HNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/resources/ebispot_harmoniser/reference/homo_sapiens-chrMT.vcf.gz \
--hm_sumstats $harm_dir/chrMT.merged_unsorted.hm \
--hm_statfile $harm_dir/chrMT.merged.log.tsv.gz $header_args \
--na_rep_in NA \
--na_rep_out NA \
--coordinate $coordinate \
--palin_mode forward;

chr=$(awk -v RS='     ' '/chromosome/{print NR; exit}' chrMT.merged_unsorted.hm)
pos=$(awk -v RS='     ' '/base_pair_location/{print NR; exit}' chrMT.merged_unsorted.hm)

head -n1 chrMT.merged_unsorted.hm > chrMT.merged.hm;
tail -n+2 chrMT.merged_unsorted.hm | sort -n -k$chr -k$pos -T$PWD >> chrMT.merged.hm
