# For different input file, two custom scripts to caculate BLUP
## 1. If the environmental phenotype files are in spearate files
add_second_column.r 1_Panicle_Architecture_2013HN.csv env 2013HN tmp_1_Panicle_Architecture_2013HN.csv
add_second_column.r 1_Panicle_Architecture_2014WH.csv env 2014WH tmp_1_Panicle_Architecture_2014WH.csv
merge_tsv.sh 1_Panicle_Architecture.tsv tmp_1_Panicle_Architecture_2013HN.csv tmp_1_Panicle_Architecture_2014WH.csv
blup.r 1_Panicle_Architecture.tsv

## 2. If a single file has phenotypes in multiple environments
blup_single_phe_multi_envs.r 3_Yield_Traits.csv BLUP/ Yield_Traits
