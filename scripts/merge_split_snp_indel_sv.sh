export SNP_VCF=$1
export INDEL_VCF=$2
export SV_VCF=$3
#bcftools sort vcf_dir/SNP.split.$SNP_VCF > tmp_vcf_dir/SNP.split.sorted.$SNP_VCF
#bcftools sort vcf_dir/INDEL.split.$INDEL_VCF > tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF
bcftools sort vcf_dir/SV.split.$SV_VCF > tmp_vcf_dir/SV.split.sorted.$SV_VCF
#bgzip -c -@ 50 tmp_vcf_dir/SNP.split.sorted.$SNP_VCF > tmp_vcf_dir/SNP.split.sorted.$SNP_VCF.gz
#bgzip -c -@ 50 tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF > tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF.gz
bgzip -c -@ 50 tmp_vcf_dir/SV.split.sorted.$SV_VCF > tmp_vcf_dir/SV.split.sorted.$SV_VCF.gz
#bcftools index -f --tbi --threads 50 tmp_vcf_dir/SNP.split.sorted.$SNP_VCF.gz
#bcftools index -f --tbi --threads 50 tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF.gz
bcftools index -f --tbi --threads 50 tmp_vcf_dir/SV.split.sorted.$SV_VCF.gz
bcftools concat -a tmp_vcf_dir/SNP.split.sorted.$SNP_VCF.gz tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF.gz tmp_vcf_dir/SV.split.sorted.$SV_VCF.gz -o vcf_dir/SNP_INDEL_SV.split.$SNP_VCF.$INDEL_VCF 
