export vcf=$1 


mkdir -p tmp_vcf_dir
cd tmp_vcf_dir
# gzip -d -k -f ../$vcf.gz
echo -e "FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr1\t1
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr2\t2
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr3\t3
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr4\t4
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr5\t5
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr6\t6
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr7\t7
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr8\t8
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr9\t9
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr10\t10
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr11\t11
FaHC_P8_Nipponbare_TEJ#0#P8.Nipponbare.TEJ.Chr12\t12" > INDEL_chr_map.txt
    bcftools annotate --rename-chrs INDEL_chr_map.txt --threads 20 -o renamed_$vcf.gz ../raw_vcf_dir/$vcf
bcftools norm -m -both renamed_$vcf.gz  -o norm_renamed_$vcf
# assign ID to vcf
awk 'BEGIN{OFS="\t"} 
/^#/ {print; next} 
{
if (NF > 1) {
	key = "INDEL-"$1"-"$2  
    count[key]++        # 计数同一位置的变异
    $3 = key"-"count[key]  
    print
	}
}
END {
	if ( NR > 0 && NF > 0) {
		printf ""
	}
}' norm_renamed_$vcf > INDEL.split.$vcf
perl -i -ne 'print unless eof && /^[\x00\s]*$/' INDEL.split.$vcf
cp INDEL.split.$vcf ../vcf_dir/


