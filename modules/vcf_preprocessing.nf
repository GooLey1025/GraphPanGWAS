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
    tag "SNP"
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
    
    # SPLIT version: normalize (split multi-allelic sites)
    ${params.bcftools} norm -m -both renamed_${sv_vcf}.gz -o norm_renamed_${sv_vcf}
    
    # Filter SV for SPLIT version (if filter script exists)
    if [ -f "${projectDir}/scripts/sv_filter.sh" ]; then
        bash ${projectDir}/scripts/sv_filter.sh norm_renamed_${sv_vcf} sv.split.temp.${sv_vcf}
    else
        cp norm_renamed_${sv_vcf} sv.split.temp.${sv_vcf}
    fi
    
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
    ${params.bcftools} view renamed_${sv_vcf}.gz > renamed_${sv_vcf}
    
    # Filter SV for UNSPLIT version (if filter script exists)
    if [ -f "${projectDir}/scripts/sv_filter.sh" ]; then
        bash ${projectDir}/scripts/sv_filter.sh renamed_${sv_vcf} sv.unsplit.temp.${sv_vcf}
    else
        cp renamed_${sv_vcf} sv.unsplit.temp.${sv_vcf}
    fi
    
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

process MERGE_VCF_FILES_SPLIT {
    tag "${analysis_type}_split"
    label 'process_very_high'
    publishDir "${params.vcf_dir}", mode: 'copy'
    
    input:
    path vcf_files
    val analysis_type
    
    output:
    path "merged.${analysis_type}_split.vcf", emit: merged_vcf
    
    script:
    def vcf_list = vcf_files.collect{ it.toString() }.join(' ')
    """
    # Sort each VCF file
    for vcf in ${vcf_list}; do
        ${params.bcftools} sort "\${vcf}" > "\${vcf}.sorted"
        ${params.bgzip} -c -@ ${task.cpus} "\${vcf}.sorted" > "\${vcf}.sorted.gz"
        ${params.bcftools} index -f --tbi --threads ${task.cpus} "\${vcf}.sorted.gz"
    done
    
    # Create list of sorted VCF files
    ls *.sorted.gz > vcf_list.txt
    
    # Merge VCF files
    ${params.bcftools} concat -a -f vcf_list.txt -o merged.${analysis_type}_split.vcf
    """
}

process MERGE_VCF_FILES_UNSPLIT {
    tag "${analysis_type}_unsplit"
    label 'process_very_high'
    publishDir "${params.vcf_dir}", mode: 'copy'
    
    input:
    path vcf_files
    val analysis_type
    
    output:
    path "merged.${analysis_type}_unsplit.vcf", emit: merged_vcf
    
    script:
    def vcf_list = vcf_files.collect{ it.toString() }.join(' ')
    """
    # Sort each VCF file
    for vcf in ${vcf_list}; do
        ${params.bcftools} sort "\${vcf}" > "\${vcf}.sorted"
        ${params.bgzip} -c -@ ${task.cpus} "\${vcf}.sorted" > "\${vcf}.sorted.gz"
        ${params.bcftools} index -f --tbi --threads ${task.cpus} "\${vcf}.sorted.gz"
    done
    
    # Create list of sorted VCF files
    ls *.sorted.gz > vcf_list.txt
    
    # Merge VCF files
    ${params.bcftools} concat -a -f vcf_list.txt -o merged.${analysis_type}_unsplit.vcf
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
    analysis_type
    
    main:
    // Process each VCF type if provided (generates both split and unsplit)
    processed_split_vcfs = Channel.empty()
    processed_unsplit_vcfs = Channel.empty()
    
    if (params.snp_vcf) {
        snp_processed = PREPARE_SNP(snp_vcf_ch)
        processed_split_vcfs = processed_split_vcfs.mix(snp_processed.split_vcf)
        processed_unsplit_vcfs = processed_unsplit_vcfs.mix(snp_processed.unsplit_vcf)
    }
    
    if (params.indel_vcf) {
        indel_processed = PREPARE_INDEL(indel_vcf_ch)
        processed_split_vcfs = processed_split_vcfs.mix(indel_processed.split_vcf)
        processed_unsplit_vcfs = processed_unsplit_vcfs.mix(indel_processed.unsplit_vcf)
    }
    
    if (params.sv_vcf) {
        sv_processed = PREPARE_SV(sv_vcf_ch)
        processed_split_vcfs = processed_split_vcfs.mix(sv_processed.split_vcf)
        processed_unsplit_vcfs = processed_unsplit_vcfs.mix(sv_processed.unsplit_vcf)
    }
    
    // Merge both split and unsplit VCF files (always generate both)
    merged_split_vcf = MERGE_VCF_FILES_SPLIT(
        processed_split_vcfs.collect(),
        analysis_type
    )
    
    merged_unsplit_vcf = MERGE_VCF_FILES_UNSPLIT(
        processed_unsplit_vcfs.collect(),
        analysis_type
    )
    
    emit:
    split_vcf = merged_split_vcf.merged_vcf
    unsplit_vcf = merged_unsplit_vcf.merged_vcf
}

