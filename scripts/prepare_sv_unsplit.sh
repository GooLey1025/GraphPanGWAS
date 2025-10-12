export vcf=$1 # Pangenie pipeline output. $vcf.gz is required to move to current dir. Note that it is compressed though environmental variable vcf is decompressed.


mkdir -p tmp_vcf_dir
cd tmp_vcf_dir
# gzip -d -k -f ../$vcf.gz
echo -e "P8.Nipponbare.TEJ.Chr1\t1
P8.Nipponbare.TEJ.Chr2\t2
P8.Nipponbare.TEJ.Chr3\t3
P8.Nipponbare.TEJ.Chr4\t4
P8.Nipponbare.TEJ.Chr5\t5
P8.Nipponbare.TEJ.Chr6\t6
P8.Nipponbare.TEJ.Chr7\t7
P8.Nipponbare.TEJ.Chr8\t8
P8.Nipponbare.TEJ.Chr9\t9
P8.Nipponbare.TEJ.Chr10\t10
P8.Nipponbare.TEJ.Chr11\t11
P8.Nipponbare.TEJ.Chr12\t12" > sv_chr_map.txt
    bcftools annotate --rename-chrs sv_chr_map.txt --threads 10 -Oz -o renamed_$vcf ../raw_vcf_dir/$vcf
# filter_SV
bash ../scripts/sv_filter.sh renamed_$vcf sv.$vcf
# assign ID to vcf
awk 'BEGIN{OFS="\t"} 
/^#/ {print; next} 
{
if (NF > 1) {
    key = "SV-"$1"-"$2  # 生成前缀 "SV-CHROM-POS"
    count[key]++        # 计数同一位置的变异
    $3 = key"-"count[key]  # 生成 "SV-CHROM-POS-COUNT" 格式的 ID
    print
	}
}
END {
	if (NR > 0 && NF > 0) {
		printf ""
	}
}' sv.$vcf > SV.unsplit.$vcf
perl -i -ne 'print unless eof && /^[\x00\s]*$/' SV.unsplit.$vcf
cp SV.unsplit.$vcf ../vcf_dir/


