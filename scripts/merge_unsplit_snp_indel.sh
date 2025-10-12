export SNP_VCF=$1
export INDEL_VCF=$2
bcftools sort vcf_dir/SNP.unsplit.$SNP_VCF > tmp_vcf_dir/SNP.unsplit.sorted.$SNP_VCF
bcftools sort vcf_dir/INDEL.unsplit.$INDEL_VCF > tmp_vcf_dir/INDEL.unsplit.sorted.$INDEL_VCF
bgzip -c -@ 50 tmp_vcf_dir/SNP.unsplit.sorted.$SNP_VCF > tmp_vcf_dir/SNP.unsplit.sorted.$SNP_VCF.gz
bgzip -c -@ 50 tmp_vcf_dir/INDEL.unsplit.sorted.$INDEL_VCF > tmp_vcf_dir/INDEL.unsplit.sorted.$INDEL_VCF.gz
bcftools index -f --tbi --threads 50 tmp_vcf_dir/SNP.unsplit.sorted.$SNP_VCF.gz
bcftools index -f --tbi --threads 50 tmp_vcf_dir/INDEL.unsplit.sorted.$INDEL_VCF.gz
bcftools concat -a tmp_vcf_dir/SNP.unsplit.sorted.$SNP_VCF.gz tmp_vcf_dir/INDEL.unsplit.sorted.$INDEL_VCF.gz -o vcf_dir/SNP_INDEL.unsplit.$SNP_VCF.$INDEL_VCF 
