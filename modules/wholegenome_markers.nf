/*
========================================================================================
    Whole Genome Markers Module
========================================================================================
    Generates LD-pruned markers for whole genome analysis
    - Converts VCF to PLINK format
    - Performs LD pruning to remove linked markers
    - Outputs pruned marker set for downstream analysis
========================================================================================
*/

workflow WHOLE_GENOME_MARKERS {
    take:
    vcf_files  // Channel of VCF files
    population_list  // Population list file
    
    main:
    // Combine VCF files with population list for each process
    vcf_with_list = vcf_files.combine(population_list)
    
    // Process each VCF file to generate LD-pruned markers
    pruned_markers = PLINK_LD_PRUNING(vcf_with_list)
    
    emit:
    pruned_markers = pruned_markers.pruned_files
}

process PLINK_LD_PRUNING {
    tag "${vcf.baseName}_LD_pruning"
    label 'process_medium'
    
    input:
    tuple path(vcf), path(population_list)
    
    output:
    tuple val("${vcf.baseName}"), 
          path("${vcf.baseName}_pruned.bed"), 
          path("${vcf.baseName}_pruned.bim"), 
          path("${vcf.baseName}_pruned.fam"),
          path("${vcf.baseName}_pruned.prune.in"),
          path("${vcf.baseName}_pruned.prune.out"),
          path("${vcf.baseName}_pruned.reheader.vcf.gz"),
          emit: pruned_files

    publishDir "${params.output_prefix}/whole_genome_markers", 
               mode: 'copy', 
               pattern: "*_pruned.*"

    script:
    """
    # Convert VCF to PLINK binary format
    ${params.plink} \\
        --vcf ${vcf} \\
        --memory ${params.plink_memory} \\
        --make-bed \\
        --out ${vcf.baseName} \\
        --allow-extra-chr \\
        --allow-no-sex
    
    # Perform LD pruning
    ${params.plink} \\
        --bfile ${vcf.baseName} \\
        --indep-pairwise 750 5 0.2 \\
        --out ${vcf.baseName}_pruned \\
        --allow-extra-chr \\
        --allow-no-sex
    
    # Extract pruned markers
    ${params.plink} \\
        --bfile ${vcf.baseName} \\
        --extract ${vcf.baseName}_pruned.prune.in \\
        --make-bed \\
        --out ${vcf.baseName}_pruned \\
        --allow-extra-chr \\
        --allow-no-sex
    
    # Extract pruned VCF using bcftools (uncompressed for reheader)
    ${params.bcftools} view \\
        -i 'ID=@${vcf.baseName}_pruned.prune.in' \\
        ${vcf} \\
        -o ${vcf.baseName}_pruned.vcf

    # Reheader with population list and compress
    ${params.bcftools} reheader -s ${population_list} -o ${vcf.baseName}_pruned.reheader.vcf ${vcf.baseName}_pruned.vcf
    ${params.bgzip} -f ${vcf.baseName}_pruned.reheader.vcf

    """
}