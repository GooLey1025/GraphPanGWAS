export SNP_VCF=$1
export INDEL_VCF=$2
bcftools sort vcf_dir/SNP.split.$SNP_VCF > tmp_vcf_dir/SNP.split.sorted.$SNP_VCF
bcftools sort vcf_dir/INDEL.split.$INDEL_VCF > tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF
bgzip -c -@ 50 tmp_vcf_dir/SNP.split.sorted.$SNP_VCF > tmp_vcf_dir/SNP.split.sorted.$SNP_VCF.gz
bgzip -c -@ 50 tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF > tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF.gz
bcftools index -f --tbi --threads 50 tmp_vcf_dir/SNP.split.sorted.$SNP_VCF.gz
bcftools index -f --tbi --threads 50 tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF.gz
bcftools concat -a tmp_vcf_dir/SNP.split.sorted.$SNP_VCF.gz tmp_vcf_dir/INDEL.split.sorted.$INDEL_VCF.gz -o vcf_dir/SNP_INDEL.split.$SNP_VCF.$INDEL_VCF 
