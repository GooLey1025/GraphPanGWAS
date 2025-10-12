/*
========================================================================================
    GWAS Analysis Module
========================================================================================
    Performs genome-wide association analysis using:
    - PLINK: VCF conversion, filtering, clumping
    - GEMMA: Kinship matrix calculation and LMM association
    - GCTA: Heritability estimation
========================================================================================
*/

process PREPROCESS_PHENOTYPE_SPLIT {
    tag "${phenotype.baseName}_split"
    label 'process_single'
    
    input:
    tuple path(vcf), path(phenotype), path(population_list)
    
    output:
    tuple path(vcf), path("processed_${phenotype.name}"), path(population_list), emit: preprocessed
    
    script:
    """
    ${params.python3} ${projectDir}/scripts/process_phenotype.py \\
        ${population_list} \\
        ${phenotype} \\
        sorted_${phenotype.name}
    
    # Add first column as FID
    ${params.awk} '{printf "%s", \$1; for(i=1; i<=NF; i++) printf "\\t%s", \$i; print "";}' \\
        sorted_${phenotype.name} > tmp_${phenotype.name}.pre_plink.csv
    
    # Format header with FID and IID
    ${params.awk} 'BEGIN {FS="\\t"; OFS="\\t"} NR==1 {\$1="FID"; \$2="IID"} {print}' \\
        tmp_${phenotype.name}.pre_plink.csv > processed_${phenotype.name}
    """
}

process PREPROCESS_PHENOTYPE_UNSPLIT {
    tag "${phenotype.baseName}_unsplit"
    label 'process_single'
    
    input:
    tuple path(vcf), path(phenotype), path(population_list)
    
    output:
    tuple path(vcf), path("processed_${phenotype.name}"), path(population_list), emit: preprocessed
    
    script:
    """
    ${params.python3} ${projectDir}/scripts/process_phenotype.py \\
        ${population_list} \\
        ${phenotype} \\
        sorted_${phenotype.name}
    
    # Add first column as FID
    ${params.awk} '{printf "%s", \$1; for(i=1; i<=NF; i++) printf "\\t%s", \$i; print "";}' \\
        sorted_${phenotype.name} > tmp_${phenotype.name}.pre_plink.csv
    
    # Format header with FID and IID
    ${params.awk} 'BEGIN {FS="\\t"; OFS="\\t"} NR==1 {\$1="FID"; \$2="IID"} {print}' \\
        tmp_${phenotype.name}.pre_plink.csv > processed_${phenotype.name}
    """
}

process PLINK_VCF_CONVERSION {
    tag "${phenotype.baseName}"
    label 'process_medium'
    
    input:
    tuple path(vcf), path(phenotype), path(population_list)
    
    output:
    tuple val("${phenotype.baseName}"), 
          path("${phenotype.baseName}.bed"), 
          path("${phenotype.baseName}.bim"), 
          path("${phenotype.baseName}.fam"),
          path(phenotype), emit: plink_files
    
    script:
    """
    ${params.plink} \\
        --vcf ${vcf} \\
        --memory ${params.plink_memory} \\
        --pheno ${phenotype} \\
        --mpheno 1 \\
        --make-bed \\
        --out ${phenotype.baseName} \\
        --allow-extra-chr \\
        --allow-no-sex
    """
}

process GEMMA_KINSHIP {
    tag "${phenotype_name}"
    label 'process_medium'
    publishDir "${params.output_prefix}/${params.analysis_type}/tmp/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.cXX.txt"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(phenotype)
    
    output:
    tuple val(phenotype_name), 
          path(bed), 
          path(bim), 
          path(fam),
          path("kinship.cXX.txt"),
          path(phenotype), emit: kinship
    
    script:
    """
    ${params.gemma} \\
        -bfile ${phenotype_name} \\
        -gk 1 \\
        -outdir . \\
        -o kinship \\
        -miss ${params.miss_threshold} \\
        -maf ${params.maf_threshold} \\
        -hwe ${params.hwe_threshold}
    """
}

process GEMMA_LMM_ASSOCIATION {
    tag "${phenotype_name}"
    label 'process_medium'
    publishDir "${params.output_prefix}/${params.analysis_type}/results/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.assoc.txt"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(kinship), path(phenotype)
    
    output:
    tuple val(phenotype_name), 
          path(bed), 
          path(bim), 
          path(fam),
          path("gemma_lmm.assoc.txt"),
          path(phenotype), emit: association
    
    script:
    """
    ${params.gemma} \\
        -bfile ${phenotype_name} \\
        -k ${kinship} \\
        -outdir . \\
        -o gemma_lmm \\
        -lmm ${params.gemma_lmm_type} \\
        -miss ${params.miss_threshold} \\
        -maf ${params.maf_threshold}
    """
}

process PLINK_CLUMPING {
    tag "${phenotype_name}"
    label 'process_low'
    publishDir "${params.output_prefix}/${params.analysis_type}/lead_markers", 
               mode: 'copy', 
               pattern: "*.clumped"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype)
    
    output:
    tuple val(phenotype_name), path("${phenotype_name}.clumped"), emit: clumped optional true
    
    script:
    """
    ${params.plink} \\
        --bfile ${phenotype_name} \\
        --clump ${assoc} \\
        --clump-p1 ${params.clump_p1} \\
        --clump-p2 ${params.clump_p2} \\
        --clump-r2 ${params.clump_r2} \\
        --clump-kb ${params.clump_kb} \\
        --allow-no-sex \\
        --allow-extra-chr \\
        --out ${phenotype_name} \\
        --clump-snp-field rs \\
        --clump-field p_wald || touch ${phenotype_name}.clumped
    """
}

process GCTA_HERITABILITY {
    tag "${phenotype_name}"
    label 'process_high'
    publishDir "${params.output_prefix}/${params.analysis_type}/results/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.hsq"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype)
    
    output:
    tuple val(phenotype_name), path("${phenotype_name}.hsq"), emit: heritability
    
    script:
    """
    # Create GRM
    ${params.gcta} \\
        --bfile ${phenotype_name} \\
        --autosome \\
        --maf ${params.maf_threshold} \\
        --make-grm \\
        --out ${phenotype_name} \\
        --thread-num ${params.gcta_threads}
    
    # Estimate heritability
    ${params.gcta} \\
        --grm ${phenotype_name} \\
        --pheno ${phenotype} \\
        --reml \\
        --out ${phenotype_name} \\
        --thread-num ${params.gcta_threads}
    """
}

/*
========================================================================================
    Main GWAS Analysis Workflow
========================================================================================
*/

workflow GWAS_ANALYSIS_SPLIT {
    take:
    gwas_input  // tuple(vcf, phenotype, population_list)
    
    main:
    // Preprocess phenotype files
    preprocessed = PREPROCESS_PHENOTYPE_SPLIT(gwas_input)
    
    // Convert VCF to PLINK binary format
    plink_files = PLINK_VCF_CONVERSION(preprocessed.preprocessed)
    
    // Calculate kinship matrix with GEMMA
    kinship = GEMMA_KINSHIP(plink_files.plink_files)
    
    // Run GEMMA LMM association
    association = GEMMA_LMM_ASSOCIATION(kinship.kinship)
    
    // Perform LD clumping with PLINK
    clumped = PLINK_CLUMPING(association.association)
    
    // Estimate heritability with GCTA
    heritability = GCTA_HERITABILITY(association.association)
    
    emit:
    assoc_files = association.association
    clumped_files = clumped.clumped
    hsq_files = heritability.heritability
}

workflow GWAS_ANALYSIS_UNSPLIT {
    take:
    gwas_input  // tuple(vcf, phenotype, population_list)
    
    main:
    // Preprocess phenotype files
    preprocessed = PREPROCESS_PHENOTYPE_UNSPLIT(gwas_input)
    
    // Convert VCF to PLINK binary format
    plink_files = PLINK_VCF_CONVERSION(preprocessed.preprocessed)
    
    // Calculate kinship matrix with GEMMA
    kinship = GEMMA_KINSHIP(plink_files.plink_files)
    
    // Run GEMMA LMM association
    association = GEMMA_LMM_ASSOCIATION(kinship.kinship)
    
    // Perform LD clumping with PLINK
    clumped = PLINK_CLUMPING(association.association)
    
    // Estimate heritability with GCTA
    heritability = GCTA_HERITABILITY(association.association)
    
    emit:
    assoc_files = association.association
    clumped_files = clumped.clumped
    hsq_files = heritability.heritability
}

