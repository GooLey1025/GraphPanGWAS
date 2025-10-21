/*
========================================================================================
    VCF Preprocessing Module
========================================================================================
    Processes SNP, INDEL, and SV VCF files:
    - Renames chromosomes
    - Normalizes variants (splits multi-allelic sites)
    - Assigns unique IDs
    - Merges VCF files based on analysis type
========================================================================================
*/

process PREPARE_SNP {
    tag "PREPARE_SNP"
    label 'process_high'
    publishDir "${params.vcf_dir}", mode: 'copy', pattern: "SNP.*.vcf"
    
    input:
    path snp_vcf
    
    output:
    path "SNP.split.${snp_vcf}", emit: split_vcf
    path "SNP.unsplit.${snp_vcf}", emit: unsplit_vcf
    
    script:
    """
    # Create chromosome mapping (shared step)
    cat > snp_chr_map.txt <<EOF
${params.chr_prefix}1\t1
${params.chr_prefix}2\t2
${params.chr_prefix}3\t3
${params.chr_prefix}4\t4
${params.chr_prefix}5\t5
${params.chr_prefix}6\t6
${params.chr_prefix}7\t7
${params.chr_prefix}8\t8
${params.chr_prefix}9\t9
${params.chr_prefix}10\t10
${params.chr_prefix}11\t11
${params.chr_prefix}12\t12
EOF
    
    # Rename chromosomes (shared step)
    ${params.bcftools} annotate --rename-chrs snp_chr_map.txt --threads ${task.cpus} \\
        -o renamed_${snp_vcf} ${snp_vcf}
    
    # SPLIT version: normalize (split multi-allelic sites)
    ${params.bcftools} norm -m -both renamed_${snp_vcf} -o norm_renamed_${snp_vcf}
    
    # Assign IDs for SPLIT version
    ${params.awk} 'BEGIN{OFS="\\t"} 
    /^#/ {print; next} 
    {
        if (NF > 1) {
            key = "SNP-"\$1"-"\$2
            count[key]++
            \$3 = key"-"count[key]
            print
        }
    }
    END {
        if (NR > 0 && NF > 0) {
            printf ""
        }
    }' norm_renamed_${snp_vcf} > SNP.split.${snp_vcf}
    
    # Assign IDs for UNSPLIT version (no normalization)
    ${params.awk} 'BEGIN{OFS="\\t"} 
    /^#/ {print; next} 
    {
        if (NF > 1) {
            key = "SNP-"\$1"-"\$2
            count[key]++
            \$3 = key"-"count[key]
            printf "%s", \$0
            printf "\\n"
        }
    }' renamed_${snp_vcf} > SNP.unsplit.${snp_vcf}
    
    # Remove trailing blank lines
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' SNP.split.${snp_vcf}
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' SNP.unsplit.${snp_vcf}
    """
}

process PREPARE_INDEL {
    tag "INDEL"
    label 'process_high'
    publishDir "${params.vcf_dir}", mode: 'copy', pattern: "INDEL.*.vcf"
    
    input:
    path indel_vcf
    
    output:
    path "INDEL.split.${indel_vcf}", emit: split_vcf
    path "INDEL.unsplit.${indel_vcf}", emit: unsplit_vcf
    
    script:
    """
    # Create chromosome mapping (shared step)
    cat > indel_chr_map.txt <<EOF
${params.chr_prefix}1\t1
${params.chr_prefix}2\t2
${params.chr_prefix}3\t3
${params.chr_prefix}4\t4
${params.chr_prefix}5\t5
${params.chr_prefix}6\t6
${params.chr_prefix}7\t7
${params.chr_prefix}8\t8
${params.chr_prefix}9\t9
${params.chr_prefix}10\t10
${params.chr_prefix}11\t11
${params.chr_prefix}12\t12
EOF
    
    # Rename chromosomes (shared step)
    ${params.bcftools} annotate --rename-chrs indel_chr_map.txt --threads ${task.cpus} \\
        -o renamed_${indel_vcf} ${indel_vcf}
    
    # SPLIT version: normalize (split multi-allelic sites)
    ${params.bcftools} norm -m -both renamed_${indel_vcf} -o norm_renamed_${indel_vcf}
    
    # Assign IDs for SPLIT version
    ${params.awk} 'BEGIN{OFS="\\t"} 
    /^#/ {print; next} 
    {
        if (NF > 1) {
            key = "INDEL-"\$1"-"\$2
            count[key]++
            \$3 = key"-"count[key]
            print
        }
    }
    END {
        if (NR > 0 && NF > 0) {
            printf ""
        }
    }' norm_renamed_${indel_vcf} > INDEL.split.${indel_vcf}
    
    # Assign IDs for UNSPLIT version (no normalization)
    ${params.awk} 'BEGIN{OFS="\\t"} 
    /^#/ {print; next} 
    {
        if (NF > 1) {
            key = "INDEL-"\$1"-"\$2
            count[key]++
            \$3 = key"-"count[key]
            printf "%s", \$0
            printf "\\n"
        }
    }' renamed_${indel_vcf} > INDEL.unsplit.${indel_vcf}
    
    # Remove trailing blank lines
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' INDEL.split.${indel_vcf}
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' INDEL.unsplit.${indel_vcf}
    """
}

process PREPARE_SV {
    tag "SV"
    label 'process_medium'
    publishDir "${params.vcf_dir}", mode: 'copy', pattern: "SV.*.vcf"
    
    input:
    path sv_vcf
    path ref_fa
    
    output:
    path "SV.split.${sv_vcf}", emit: split_vcf
    path "SV.unsplit.${sv_vcf}", emit: unsplit_vcf
    
    script:
    """
    # Create chromosome mapping (SV uses different prefix - shared step)
    cat > sv_chr_map.txt <<EOF
P8.Nipponbare.TEJ.Chr1\t1
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
P8.Nipponbare.TEJ.Chr12\t12
EOF
    
    # Rename chromosomes (shared step)
    ${params.bcftools} annotate --rename-chrs sv_chr_map.txt --threads ${task.cpus} \\
        -Oz -o renamed_${sv_vcf}.gz ${sv_vcf}

    # left-align the vcf
    ${params.bcftools} norm -f ${ref_fa} -Oz -o LeftAligned_renamed_${sv_vcf}.gz renamed_${sv_vcf}.gz

    # SPLIT version: normalize (split multi-allelic sites)
    ${params.bcftools} norm -m -both LeftAligned_renamed_${sv_vcf}.gz -o norm_LeftAligned_renamed_${sv_vcf}
    
    # Filter SV for SPLIT version 
    bash ${projectDir}/scripts/sv_filter.sh norm_LeftAligned_renamed_${sv_vcf} sv.split.temp.${sv_vcf}
    
    # Assign IDs for SPLIT version
    ${params.awk} 'BEGIN{OFS="\\t"} 
    /^#/ {print; next} 
    {
        if (NF > 1) {
            key = "SV-"\$1"-"\$2
            count[key]++
            \$3 = key"-"count[key]
            print
        }
    }
    END {
        if (NR > 0 && NF > 0) {
            printf ""
        }
    }' sv.split.temp.${sv_vcf} > SV.split.${sv_vcf}
    
    # UNSPLIT version: extract from gzipped renamed file
    ${params.bcftools} view LeftAligned_renamed_${sv_vcf}.gz > LeftAligned_renamed_${sv_vcf}
    
    # Filter SV for UNSPLIT version 
    bash ${projectDir}/scripts/sv_filter.sh LeftAligned_renamed_${sv_vcf} sv.unsplit.temp.${sv_vcf}

    
    # Assign IDs for UNSPLIT version (no normalization)
    ${params.awk} 'BEGIN{OFS="\\t"} 
    /^#/ {print; next} 
    {
        if (NF > 1) {
            key = "SV-"\$1"-"\$2
            count[key]++
            \$3 = key"-"count[key]
            printf "%s", \$0
            printf "\\n"
        }
    }' sv.unsplit.temp.${sv_vcf} > SV.unsplit.${sv_vcf}
    
    # Remove trailing blank lines
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' SV.split.${sv_vcf}
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' SV.unsplit.${sv_vcf}
    """
}


process MERGE_SNP_INDEL_SPLIT {
    tag "SNP_INDEL_split"
    label 'process_very_high'
    publishDir "${params.vcf_dir}", mode: 'copy'
    
    input:
    path snp_vcf
    path indel_vcf
    
    output:
    path "SNP_INDEL.split.*.vcf", emit: merged_vcf
    
    script:
    def output_parts = []
    output_parts.add(snp_vcf.name.replaceFirst(/^SNP\.split\./, ''))
    output_parts.add(indel_vcf.name.replaceFirst(/^INDEL\.split\./, ''))
    def output_name = "SNP_INDEL.split." + output_parts.join('.')
    
    """
    # Create tmp directory for intermediate files
    mkdir -p tmp_vcf_dir
    
    # Process SNP file
    echo "Processing SNP file: ${snp_vcf}"
    ${params.bcftools} sort "${snp_vcf}" -o "tmp_vcf_dir/SNP.split.pos_sorted.\$(basename ${snp_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/SNP.split.pos_sorted.\$(basename ${snp_vcf})" | sort > "tmp_vcf_dir/snp_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/snp_samples.txt" "tmp_vcf_dir/SNP.split.pos_sorted.\$(basename ${snp_vcf})" -o "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf})" > "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf}).gz"
    
    # Process INDEL file
    echo "Processing INDEL file: ${indel_vcf}"
    ${params.bcftools} sort "${indel_vcf}" -o "tmp_vcf_dir/INDEL.split.pos_sorted.\$(basename ${indel_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/INDEL.split.pos_sorted.\$(basename ${indel_vcf})" | sort > "tmp_vcf_dir/indel_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/indel_samples.txt" "tmp_vcf_dir/INDEL.split.pos_sorted.\$(basename ${indel_vcf})" -o "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf})" > "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf}).gz"
    
    # Merge files
    ${params.bcftools} concat -a tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf}).gz tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf}).gz -o ${output_name}
    """
}

process MERGE_SNP_INDEL_SV_SPLIT {
    tag "SNP_INDEL_SV_split"
    label 'process_very_high'
    publishDir "${params.vcf_dir}", mode: 'copy'
    
    input:
    path snp_vcf
    path indel_vcf
    path sv_vcf
    
    output:
    path "SNP_INDEL_SV.split.*.vcf", emit: merged_vcf
    
    script:
    def output_parts = []
    output_parts.add(snp_vcf.name.replaceFirst(/^SNP\.split\./, ''))
    output_parts.add(indel_vcf.name.replaceFirst(/^INDEL\.split\./, ''))
    output_parts.add(sv_vcf.name.replaceFirst(/^SV\.split\./, ''))
    def output_name = "SNP_INDEL_SV.split." + output_parts.join('.')
    
    """
    # Create tmp directory for intermediate files
    mkdir -p tmp_vcf_dir
    
    # Process SNP file
    echo "Processing SNP file: ${snp_vcf}"
    ${params.bcftools} sort "${snp_vcf}" -o "tmp_vcf_dir/SNP.split.pos_sorted.\$(basename ${snp_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/SNP.split.pos_sorted.\$(basename ${snp_vcf})" | sort > "tmp_vcf_dir/snp_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/snp_samples.txt" "tmp_vcf_dir/SNP.split.pos_sorted.\$(basename ${snp_vcf})" -o "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf})" > "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf}).gz"
    
    # Process INDEL file
    echo "Processing INDEL file: ${indel_vcf}"
    ${params.bcftools} sort "${indel_vcf}" -o "tmp_vcf_dir/INDEL.split.pos_sorted.\$(basename ${indel_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/INDEL.split.pos_sorted.\$(basename ${indel_vcf})" | sort > "tmp_vcf_dir/indel_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/indel_samples.txt" "tmp_vcf_dir/INDEL.split.pos_sorted.\$(basename ${indel_vcf})" -o "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf})" > "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf}).gz"
    
    # Process SV file
    echo "Processing SV file: ${sv_vcf}"
    ${params.bcftools} sort "${sv_vcf}" -o "tmp_vcf_dir/SV.split.pos_sorted.\$(basename ${sv_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/SV.split.pos_sorted.\$(basename ${sv_vcf})" | sort > "tmp_vcf_dir/sv_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/sv_samples.txt" "tmp_vcf_dir/SV.split.pos_sorted.\$(basename ${sv_vcf})" -o "tmp_vcf_dir/SV.split.sorted.\$(basename ${sv_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/SV.split.sorted.\$(basename ${sv_vcf})" > "tmp_vcf_dir/SV.split.sorted.\$(basename ${sv_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/SV.split.sorted.\$(basename ${sv_vcf}).gz"
    
    # Merge files
    ${params.bcftools} concat -a tmp_vcf_dir/SNP.split.sorted.\$(basename ${snp_vcf}).gz tmp_vcf_dir/INDEL.split.sorted.\$(basename ${indel_vcf}).gz tmp_vcf_dir/SV.split.sorted.\$(basename ${sv_vcf}).gz -o ${output_name}
    """
}

process MERGE_SNP_INDEL_UNSPLIT {
    tag "SNP_INDEL_unsplit"
    label 'process_very_high'
    publishDir "${params.vcf_dir}", mode: 'copy'
    
    input:
    path snp_vcf
    path indel_vcf
    
    output:
    path "SNP_INDEL.unsplit.*.vcf", emit: merged_vcf
    
    script:
    def output_parts = []
    output_parts.add(snp_vcf.name.replaceFirst(/^SNP\.unsplit\./, ''))
    output_parts.add(indel_vcf.name.replaceFirst(/^INDEL\.unsplit\./, ''))
    def output_name = "SNP_INDEL.unsplit." + output_parts.join('.')
    
    """
    # Create tmp directory for intermediate files
    mkdir -p tmp_vcf_dir
    
    # Process SNP file
    echo "Processing SNP file: ${snp_vcf}"
    ${params.bcftools} sort "${snp_vcf}" -o "tmp_vcf_dir/SNP.unsplit.pos_sorted.\$(basename ${snp_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/SNP.unsplit.pos_sorted.\$(basename ${snp_vcf})" | sort > "tmp_vcf_dir/snp_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/snp_samples.txt" "tmp_vcf_dir/SNP.unsplit.pos_sorted.\$(basename ${snp_vcf})" -o "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf})" > "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf}).gz"
    
    # Process INDEL file
    echo "Processing INDEL file: ${indel_vcf}"
    ${params.bcftools} sort "${indel_vcf}" -o "tmp_vcf_dir/INDEL.unsplit.pos_sorted.\$(basename ${indel_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/INDEL.unsplit.pos_sorted.\$(basename ${indel_vcf})" | sort > "tmp_vcf_dir/indel_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/indel_samples.txt" "tmp_vcf_dir/INDEL.unsplit.pos_sorted.\$(basename ${indel_vcf})" -o "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf})" > "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf}).gz"
    
    # Merge files
    ${params.bcftools} concat -a tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf}).gz tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf}).gz -o ${output_name}
    """
}

process MERGE_SNP_INDEL_SV_UNSPLIT {
    tag "SNP_INDEL_SV_unsplit"
    label 'process_very_high'
    publishDir "${params.vcf_dir}", mode: 'copy'
    
    input:
    path snp_vcf
    path indel_vcf
    path sv_vcf
    
    output:
    path "SNP_INDEL_SV.unsplit.*.vcf", emit: merged_vcf
    
    script:
    def output_parts = []
    output_parts.add(snp_vcf.name.replaceFirst(/^SNP\.unsplit\./, ''))
    output_parts.add(indel_vcf.name.replaceFirst(/^INDEL\.unsplit\./, ''))
    output_parts.add(sv_vcf.name.replaceFirst(/^SV\.unsplit\./, ''))
    def output_name = "SNP_INDEL_SV.unsplit." + output_parts.join('.')
    
    """
    # Create tmp directory for intermediate files
    mkdir -p tmp_vcf_dir
    
    # Process SNP file
    echo "Processing SNP file: ${snp_vcf}"
    ${params.bcftools} sort "${snp_vcf}" -o "tmp_vcf_dir/SNP.unsplit.pos_sorted.\$(basename ${snp_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/SNP.unsplit.pos_sorted.\$(basename ${snp_vcf})" | sort > "tmp_vcf_dir/snp_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/snp_samples.txt" "tmp_vcf_dir/SNP.unsplit.pos_sorted.\$(basename ${snp_vcf})" -o "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf})" > "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf}).gz"
    
    # Process INDEL file
    echo "Processing INDEL file: ${indel_vcf}"
    ${params.bcftools} sort "${indel_vcf}" -o "tmp_vcf_dir/INDEL.unsplit.pos_sorted.\$(basename ${indel_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/INDEL.unsplit.pos_sorted.\$(basename ${indel_vcf})" | sort > "tmp_vcf_dir/indel_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/indel_samples.txt" "tmp_vcf_dir/INDEL.unsplit.pos_sorted.\$(basename ${indel_vcf})" -o "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf})" > "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf}).gz"
    
    # Process SV file
    echo "Processing SV file: ${sv_vcf}"
    ${params.bcftools} sort "${sv_vcf}" -o "tmp_vcf_dir/SV.unsplit.pos_sorted.\$(basename ${sv_vcf})"
    ${params.bcftools} query -l "tmp_vcf_dir/SV.unsplit.pos_sorted.\$(basename ${sv_vcf})" | sort > "tmp_vcf_dir/sv_samples.txt"
    ${params.bcftools} view -S "tmp_vcf_dir/sv_samples.txt" "tmp_vcf_dir/SV.unsplit.pos_sorted.\$(basename ${sv_vcf})" -o "tmp_vcf_dir/SV.unsplit.sorted.\$(basename ${sv_vcf})"
    ${params.bgzip} -c -@ ${task.cpus} "tmp_vcf_dir/SV.unsplit.sorted.\$(basename ${sv_vcf})" > "tmp_vcf_dir/SV.unsplit.sorted.\$(basename ${sv_vcf}).gz"
    ${params.bcftools} index -f --tbi --threads ${task.cpus} "tmp_vcf_dir/SV.unsplit.sorted.\$(basename ${sv_vcf}).gz"
    
    # Merge files
    ${params.bcftools} concat -a tmp_vcf_dir/SNP.unsplit.sorted.\$(basename ${snp_vcf}).gz tmp_vcf_dir/INDEL.unsplit.sorted.\$(basename ${indel_vcf}).gz tmp_vcf_dir/SV.unsplit.sorted.\$(basename ${sv_vcf}).gz -o ${output_name}
    """
}

/*
========================================================================================
    Main VCF Preprocessing Workflow
========================================================================================
*/

workflow VCF_PREPROCESSING {
    take:
    snp_vcf_ch
    indel_vcf_ch
    sv_vcf_ch
    ref_fa
    
    main:
    // Step 1: PREPARE - Process each VCF type if provided (generates both split and unsplit)
    snp_split_vcf = Channel.empty()
    snp_unsplit_vcf = Channel.empty()
    indel_split_vcf = Channel.empty()
    indel_unsplit_vcf = Channel.empty()
    sv_split_vcf = Channel.empty()
    sv_unsplit_vcf = Channel.empty()
    
    if (params.snp_vcf) {
        snp_processed = PREPARE_SNP(snp_vcf_ch)
        snp_split_vcf = snp_processed.split_vcf
        snp_unsplit_vcf = snp_processed.unsplit_vcf
    }
    
    if (params.indel_vcf) {
        indel_processed = PREPARE_INDEL(indel_vcf_ch)
        indel_split_vcf = indel_processed.split_vcf
        indel_unsplit_vcf = indel_processed.unsplit_vcf
    }
    
    if (params.sv_vcf) {
        sv_processed = PREPARE_SV(sv_vcf_ch, ref_fa)
        sv_split_vcf = sv_processed.split_vcf
        sv_unsplit_vcf = sv_processed.unsplit_vcf
    }
    
    // Step 2: MERGE - Generate all 6 types of VCF files for GWAS
    // Following the logic in run_pipeline.sh VCF_WAY_PAIRS
    
    // Create channels for different VCF combinations
    all_split_vcfs = Channel.empty()
    all_unsplit_vcfs = Channel.empty()
    
    // Type 1: SNP only (split and unsplit)
    // No merge needed for single SNP file - directly use PREPARE_SNP output
    if (params.snp_vcf) {
        all_split_vcfs = all_split_vcfs.concat(snp_split_vcf)
        all_unsplit_vcfs = all_unsplit_vcfs.concat(snp_unsplit_vcf)
    }
    
    // Type 2: SNP_INDEL (merge SNP + INDEL, split and unsplit)
    if (params.snp_vcf && params.indel_vcf) {
        snp_indel_split_merged = MERGE_SNP_INDEL_SPLIT(snp_split_vcf, indel_split_vcf)
        snp_indel_unsplit_merged = MERGE_SNP_INDEL_UNSPLIT(snp_unsplit_vcf, indel_unsplit_vcf)
        
        all_split_vcfs = all_split_vcfs.concat(snp_indel_split_merged.merged_vcf)
        all_unsplit_vcfs = all_unsplit_vcfs.concat(snp_indel_unsplit_merged.merged_vcf)
    }
    
    // Type 3: SNP_INDEL_SV (merge SNP + INDEL + SV, split and unsplit)
    if (params.snp_vcf && params.indel_vcf && params.sv_vcf) {
        snp_indel_sv_split_merged = MERGE_SNP_INDEL_SV_SPLIT(snp_split_vcf, indel_split_vcf, sv_split_vcf)
        snp_indel_sv_unsplit_merged = MERGE_SNP_INDEL_SV_UNSPLIT(snp_unsplit_vcf, indel_unsplit_vcf, sv_unsplit_vcf)
        
        all_split_vcfs = all_split_vcfs.concat(snp_indel_sv_split_merged.merged_vcf)
        all_unsplit_vcfs = all_unsplit_vcfs.concat(snp_indel_sv_unsplit_merged.merged_vcf)
    }
    
    emit:
    split_vcf = all_split_vcfs
    unsplit_vcf = all_unsplit_vcfs
}

