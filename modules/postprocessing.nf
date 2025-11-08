/*
========================================================================================
    Post-processing Module
========================================================================================
    Extracts and summarizes results:
    - GCTA heritability estimates
    - Lead marker SNPs from clumping
========================================================================================
*/

process EXTRACT_HERITABILITY {
    tag "heritability_${way}"
    label 'process_single'
    publishDir "${params.heritability_dir}", 
               mode: 'copy'
    
    input:
    tuple val(way), path(hsq_files)
    val output_prefix
    
    output:
    path "${output_prefix}.${way}.gcta.tsv", emit: heritability_summary
    
    script:
    """
    # Create output file with header
    echo -e "Phenotype\\tHeritability" > ${output_prefix}.${way}.gcta.tsv
    
    # Extract heritability from each .hsq file
    for file in *.hsq; do
        if [ -f "\${file}" ]; then
            # Get phenotype name (remove .hsq extension)
            phenotype=\$(basename "\${file}" .hsq)
            
            # Extract V(G)/Vp value from line 5, column 2
            if [ \$(wc -l < "\${file}") -ge 5 ]; then
                value=\$(sed -n '5p' "\${file}" | awk '{print \$2}')
                if [ ! -z "\${value}" ]; then
                    echo -e "\${phenotype}\\t\${value}" >> ${output_prefix}.${way}.gcta.tsv
                fi
            fi
        fi
    done
    
    echo "Heritability extraction for ${way} complete."
    """
}

process EXTRACT_LDAK_HERITABILITY {
    tag "ldak_heritability_${way}"
    label 'process_single'
    publishDir "${params.heritability_dir}", 
               mode: 'copy'
    
    input:
    tuple val(way), val(reml_file_paths)
    val output_prefix
    
    output:
    path "${output_prefix}.${way}.ldak.tsv", emit: ldak_heritability_summary
    
    script:
    def reml_paths = reml_file_paths instanceof List ? reml_file_paths.join(' ') : reml_file_paths
    """
    # Create output file with header
    echo -e "Phenotype\\tHeritability\\tSE\\tP_value" > ${output_prefix}.${way}.ldak.tsv
    
    # Extract heritability from each .reml file
    # reml_file_paths contains absolute paths to all reml files for this way
    for file in ${reml_paths}; do
        if [ -f "\${file}" ]; then
            # Get phenotype name (remove .reml extension, LDAK suffix, and processed_ prefix, then add .tsv)
            phenotype=\$(basename "\${file}" .reml | sed 's/_LDAK-Thin\$//' | sed 's/^processed_//')
            phenotype="\${phenotype}.tsv"
            
            # Extract Her_All value (total heritability) - works for both single GRM and MGRM
            # Find the line starting with "Her_All" and extract columns 2 (heritability) and 3 (SE)
            # Extract LRT_P from line 17, column 2 (p-value)
            her_all_line=\$(grep "^Her_All" "\${file}")
            if [ ! -z "\${her_all_line}" ]; then
                heritability=\$(echo "\${her_all_line}" | awk '{print \$2}')
                se=\$(echo "\${her_all_line}" | awk '{print \$3}')
                p_value=\$(sed -n '17p' "\${file}" | awk '{print \$2}')
                
                if [ ! -z "\${heritability}" ] && [ "\${heritability}" != "NA" ]; then
                    echo -e "\${phenotype}\\t\${heritability}\\t\${se}\\t\${p_value}" >> ${output_prefix}.${way}.ldak.tsv
                fi
            fi
        fi
    done
    
    echo "LDAK heritability extraction for ${way} complete."
    """
}

process GENERATE_HERITABILITY_TABLE {
    tag "heritability_table"
    label 'process_single'
    publishDir "${params.heritability_dir}", 
               mode: 'copy'
    
    input:
    path ldak_tsv_files  // Multiple TSV files from EXTRACT_LDAK_HERITABILITY
    val output_prefix
    val heritability_dir
    
    output:
    path "${output_prefix}.heritability_summary.xlsx", emit: heritability_table
    
    script:
    """
    # The files are already in the working directory or will be symlinked
    # Python script needs access to the output directory for REML files
    # Pass the full path to the output directory
    ${params.python3} ${projectDir}/scripts/combine_heritability_table.py \\
        ${projectDir}/${params.output_prefix} \\
        ${output_prefix} \\
        ${output_prefix}.heritability_summary.xlsx
    
    echo "Heritability table generation complete."
    """
}

process EXTRACT_LEAD_MARKERS {
    tag "lead_markers_${way}"
    label 'process_single'
    publishDir "${params.lead_marker_dir}", 
               mode: 'copy'
    
    input:
    tuple val(way), path(clumped_files)
    val output_prefix
    
    output:
    tuple val(way), path("${output_prefix}.${way}.ld_marker.tsv"), emit: lead_markers_summary
    
    script:
    """
    # Create output file with header
    echo -e "Marker_ID" > ${output_prefix}.${way}.ld_marker.tsv
    
    # Extract lead marker IDs from each .clumped file
    for file in *.clumped; do
        if [ -f "\${file}" ] && [ -s "\${file}" ]; then
            # Skip header line and extract SNP IDs (column 3)
            # Check if file has more than just header
            if [ \$(wc -l < "\${file}") -gt 1 ]; then
                sed -n '2,\$p' "\${file}" | awk 'NF > 0 {print \$3}' >> ${output_prefix}.${way}.ld_marker.tsv
            fi
        fi
    done
    
    # Remove duplicates while maintaining order
    sort -u ${output_prefix}.${way}.ld_marker.tsv > temp.tsv
    mv temp.tsv ${output_prefix}.${way}.ld_marker.tsv
    
    echo "Lead marker extraction for ${way} complete."
    """
}

process EXTRACT_LD_MARKER_VCF {
    tag "ld_vcf_${way}"
    label 'process_medium'
    publishDir "${params.output_prefix}/${way}/ld_filtered_vcf", 
               mode: 'copy'
    
    input:
    tuple val(way), path(ld_marker_tsv), path(vcf), path(population_list)
    
    output:
    tuple val(way), path("${vcf.baseName}.ld.reheader.vcf.gz"), emit: ld_filtered_vcf
    
    script:
    """
    # Extract VCF variants matching LD marker IDs (uncompressed for reheader)
    ${params.bcftools} view -i 'ID=@${ld_marker_tsv}' ${vcf} > ${vcf.baseName}.ld.vcf

    # Reheader with population list and compress
    ${params.bcftools} reheader -s ${population_list} -o ${vcf.baseName}.ld.reheader.vcf ${vcf.baseName}.ld.vcf
    ${params.bgzip} -f ${vcf.baseName}.ld.reheader.vcf

    echo "LD marker VCF extraction and reheader for ${way} complete."
    """
}

process GENERATE_MANHATTAN_PLOT {
    tag "${phenotype_name}"
    label 'process_low'
    publishDir "${params.output_prefix}/${params.analysis_type}/plots", 
               mode: 'copy',
               pattern: "*.png"
    
    input:
    tuple val(phenotype_name), path(assoc_file)
    
    output:
    path "${phenotype_name}_manhattan.png", emit: manhattan_plot optional true
    
    when:
    params.generate_plots
    
    script:
    """
    if [ -f "${projectDir}/scripts/manhattan_cmplot.R" ]; then
        Rscript ${projectDir}/scripts/manhattan_cmplot.R \\
            ${assoc_file} \\
            ${phenotype_name}_manhattan.png
    else
        echo "Manhattan plot script not found. Skipping..."
        touch ${phenotype_name}_manhattan.png
    fi
    """
}

/*
========================================================================================
    Summary Statistics Process
========================================================================================
*/

process GENERATE_SUMMARY_REPORT {
    tag "summary_report"
    label 'process_single'
    publishDir "${params.output_prefix}/${params.analysis_type}", 
               mode: 'copy'
    
    input:
    path heritability_file
    path lead_markers_file
    val output_prefix
    val analysis_type
    
    output:
    path "GWAS_summary_report.txt", emit: summary_report
    
    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    cat > GWAS_summary_report.txt <<EOF
========================================
GWAS Analysis Summary Report
========================================
Analysis Type: ${analysis_type}
Output Prefix: ${output_prefix}
Generated: \$(date)

========================================
Heritability Estimates
========================================
\$(cat ${heritability_file})

========================================
Lead Markers Summary
========================================
Total unique lead markers: \$(tail -n +2 ${lead_markers_file} | wc -l)

Top 10 lead markers:
\$(head -n 11 ${lead_markers_file})

========================================
Output Files Location
========================================
- Results: ${params.output_prefix}/${params.analysis_type}/results/
- Lead markers: ${params.output_prefix}/${params.analysis_type}/lead_markers/
- Heritability: ${params.output_prefix}/${params.analysis_type}/heritability/
- Pipeline info: ${params.output_prefix}/pipeline_info/

========================================
EOF

    echo "Summary report generated successfully."
    """
}

