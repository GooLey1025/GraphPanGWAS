export list=532rice.list # move to current directory.  The list of population samples
export SNP_VCF=532rice.snp.impute.vcf # GATK-pipeline output. Move to ./raw_vcf_dir/
export INDEL_VCF=532rice.indel.impute.vcf # GATK-pipeline output. Move to ./raw_vcf_dir/
export SV_VCF=532rice.pangenie.multi.impute.vcf # Pangenie-pipeline output. Move to ./raw_vcf_dir/
export P=532rice # custom output_dir_name.
export phenotypes_dir="./phenotypes/532rice" # prepare your phenotypes.tsv files. Note that suffix is .tsv. 
export error_log="$P.error_log.txt"
export way=SNP_INDEL_SV_split # The name of output_dir.
export INPUT_VCF=SNP_INDEL_SV.split.532rice.snp.impute.vcf.532rice.indel.impute.vcf # VCFs that have been processed by scripts below (e.g prepare_snp_split.sh) and they have been moved to ./vcf_dir
mkdir -p tmp_vcf_dir
mkdir -p raw_vcf_dir 
mkdir -p vcf_dir # There are vcfs which have been processed, can be input to gwas.sh

# Here are some examples to process vcf including rename "CHROM name" and "Assigning ID". Perhaps you could write you own script to process your VCF that depends on what you are going to do. Anyway, make sure they have unique ID, and it is recommended that #CHROM column is concise.
bash scripts/prepare_sv_split.sh $SV_VCF
bash scripts/prepare_sv_unsplit.sh $SV_VCF
bash scripts/prepare_snp_split.sh $SNP_VCF
bash scripts/prepare_snp_unsplit.sh $SNP_VCF
bash scripts/prepare_indel_unsplit.sh $INDEL_VCF
bash scripts/prepare_indel_split.sh $INDEL_VCF
bash scripts/merge_split_snp_indel.sh $SNP_VCF $INDEL_VCF
bash scripts/merge_unsplit_snp_indel.sh $SNP_VCF $INDEL_VCF
bash scripts/merge_split_snp_indel_sv.sh $SNP_VCF $INDEL_VCF $SV_VCF
bash scripts/merge_unsplit_snp_indel_sv.sh $SNP_VCF $INDEL_VCF $SV_VCF

> "$error_log"
# Make sure that the phenotype file has a header, with the first column being the ID and the rest being the values. It can be non-sequenced. Two columns. And the file name is a distinct phenotype name.
check_and_run() {
    export file=$1
    # Check if the file has at least two columns (ID and phenotype columns)
    export num_columns=$(head -n 1 "$file" | awk -F'\t' '{print NF}')
    
    if [ "$num_columns" -ne 2 ]; then
        echo "Error: $file does not have exactly two columns (ID and phenotype). Skipping..." >> "$error_log"
        return 1
    fi

    # Check if the first column is 'ID'
    first_column_name=$(head -n 1 "$file" | cut -f1)
    
    if [ "$first_column_name" != "ID" ]; then
        echo "Error: $file does not have 'ID' as the first column. Skipping..." >> "$error_log"
        return 1
    fi

    # If the format is correct, run the command
    echo "Running sv.multi.gwas.sh on $file..."
    bash gwas.sh "$file" $list $INPUT_VCF $P $way
}

export -f check_and_run  # Export function for parallel execution

find "$phenotypes_dir" -type f -name "*.tsv" |  parallel -j 5 check_and_run

if [ -s "$error_log" ]; then
    echo "Some files were skipped due to errors. Please check the error log at $error_log."
else
    echo "All files processed successfully."
fi

# extract GCTA variance explained
mkdir -p hertiability
> hertiability/$P.$way.gcta.tsv
find $P/$way/results/*.tsv -type f -name "*.hsq" | while read file; do
    basename=$(basename "$file" .tsv.hsq)
    value=$(sed -n '5p' "$file" | cut -f2)
    echo -e "$basename\t$value" >> hertiability/$P.$way.gcta.tsv
done

# extract Marker
mkdir -p LD_marker
> Lead_marker/$P.$way.ld_marker.tsv
find $P/$way/lead_markers/*.clumped -type f -name "*.clumped" | while read file; do
    values=$(sed -n '2,$p' "$file" | awk '{print $3}')
    for value in $values; do
        echo -e "$value" >> Lead_marker/$P.$way.ld_marker.tsv
    done
done


