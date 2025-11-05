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
    tag "${phenotype.baseName}_${way}_split"
    label 'process_single'
    
    input:
    tuple path(vcf), path(phenotype), path(population_list), val(way)
    
    output:
    tuple path(vcf), path("processed_${phenotype.name}"), path(population_list), val(way), emit: preprocessed
    
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
    tag "${phenotype.baseName}_${way}_unsplit"
    label 'process_single'
    
    input:
    tuple path(vcf), path(phenotype), path(population_list), val(way)
    
    output:
    tuple path(vcf), path("processed_${phenotype.name}"), path(population_list), val(way), emit: preprocessed
    
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
    tag "${phenotype.baseName}_${way}"
    label 'process_medium'
    
    input:
    tuple path(vcf), path(phenotype), path(population_list), val(way)
    
    output:
    tuple val("${phenotype.baseName}"), 
          path("${phenotype.baseName}.bed"), 
          path("${phenotype.baseName}.bim"), 
          path("${phenotype.baseName}.fam"),
          path(phenotype),
          val(way), emit: plink_files
    
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
    tag "${phenotype_name}_${way}"
    label 'process_medium'
    publishDir "${params.output_prefix}/${way}/tmp/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.cXX.txt"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(phenotype), val(way)
    
    output:
    tuple val(phenotype_name), 
          path(bed), 
          path(bim), 
          path(fam),
          path("kinship.cXX.txt"),
          path(phenotype),
          val(way), emit: kinship
    
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
    tag "${phenotype_name}_${way}"
    label 'process_medium'
    publishDir "${params.output_prefix}/${way}/results/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.assoc.txt"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(kinship), path(phenotype), val(way)
    
    output:
    tuple val(phenotype_name), 
          path(bed), 
          path(bim), 
          path(fam),
          path("gemma_lmm.assoc.txt"),
          path(phenotype),
          val(way), emit: association
    
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
    tag "${phenotype_name}_${way}"
    label 'process_low'
    publishDir "${params.output_prefix}/${way}/lead_markers", 
               mode: 'copy', 
               pattern: "*.clumped"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype), val(way)
    
    output:
    tuple val(phenotype_name), path("${phenotype_name}.clumped"), val(way), emit: clumped
    
    script:
    """
    set +e  # Disable exit on error for clumping
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
        --clump-field p_wald
    
    # If clumping failed or produced no output, create empty clumped file
    if [ ! -f ${phenotype_name}.clumped ]; then
        touch ${phenotype_name}.clumped
    fi
    set -e  # Re-enable exit on error
    """
}

process GCTA_HERITABILITY {
    tag "${phenotype_name}_${way}"
    label 'process_high'
    publishDir "${params.output_prefix}/${way}/results/${phenotype_name}", 
               mode: 'copy', 
               pattern: "*.hsq"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype), val(way)
    
    output:
    tuple val(phenotype_name), path("${phenotype_name}.hsq"), val(way), emit: heritability
    
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

process RECODE_VCF_FOR_LDAK {
    tag "RECODE_${vcf.baseName}_${way}"
    label 'process_medium'
    
    input:
    tuple path(vcf), val(way)
    
    output:
    tuple path("ldak_${vcf.name}"), val(way), emit: recoded_vcf
    
    script:
    """
    # Recode multi-character alleles (REF or ALT length > 1) to X and Y for LDAK
    # This is only needed for INDEL, SV, and merged VCFs (not SNP-only)
    # Use ldak_ prefix to clearly identify this is for LDAK and avoid overwriting input
    ${params.awk} 'BEGIN{OFS="\\t"}
    /^#/ {print; next}
    {
      ref=\$4; alt=\$5;
      if(length(ref)>1 || length(alt)>1){
        \$4="X"; \$5="Y"
      }
      print
    }' ${vcf} > ldak_${vcf.name}.tmp
    
    # Remove trailing blank lines from temporary file
    ${params.perl} -i -ne 'print unless eof && /^[\\x00\\s]*\$/' ldak_${vcf.name}.tmp
    
    # Copy to final output (keep tmp file for debugging)
    cp ldak_${vcf.name}.tmp ldak_${vcf.name}
    """
}

process LDAK_PREPARE_GRM {
    tag "LDAK_GRM_${way}"
    label 'process_high'
    publishDir "${params.output_prefix}/${way}/grm", 
               mode: 'copy', 
               pattern: "*.grm.*"
    
    input:
    tuple path(vcf), val(way)
    
    output:
    tuple val(way), 
          path("*_WG.bed"), 
          path("*_WG.bim"), 
          path("*_WG.fam"), 
          path("*_LDAK-Thin.grm.bin"), 
          path("*_LDAK-Thin.grm.id"), 
          path("*_LDAK-Thin.grm.details"),
          path("*_LDAK-Thin.grm.adjust"),
          emit: grm_files
    
    script:
    """
    # Convert VCF to PLINK binary format (without phenotype)
    plink --vcf ${vcf} --memory 10000 --make-bed --out ${vcf.baseName}_WG --allow-extra-chr --allow-no-sex
    
    # Cut weights regions
    ${params.ldak} --bfile ${vcf.baseName}_WG --cut-weights cutweights --window-prune 0.98
    
    # Calculate weights
    ${params.ldak} --bfile ${vcf.baseName}_WG --calc-weights-all cutweights --max-threads ${params.ldak_threads} > weights.out
    
    # Calculate GRM
    ${params.ldak} --calc-kins-direct ${vcf.baseName}_WG_LDAK-Thin --bfile ${vcf.baseName}_WG --weights cutweights/weights.all --power -.5
    """
}

process LDAK_CALCULATE_HERITABILITY {
    tag "LDAK_${phenotype.baseName}_${way}"
    label 'process_medium'
    publishDir "${params.output_prefix}/${way}/results/${phenotype.baseName}", 
               mode: 'copy', 
               pattern: "*.reml"
    
    input:
    tuple path(vcf), path(phenotype), path(population_list), val(way), val(way2), path(bed), path(bim), path(fam), path(grm_bin), path(grm_id), path(grm_details), path(grm_adjust)
    
    output:
    tuple val("${phenotype.baseName}"), path("*.reml"), val(way), emit: heritability
    
    script:
    """
    # Extract phenotype name from filename
    phenotype_name=\$(basename ${phenotype} | sed 's/processed_//')
    
    # Get GRM prefix (remove .grm.bin extension)
    # Use the full path to ensure GRM files are found in the same directory
    grm_prefix=\$(dirname ${grm_bin})/\$(basename ${grm_bin} .grm.bin)
    
    # Calculate heritability using pre-computed GRM
    # GRM details file is automatically found by LDAK using the prefix
    ${params.ldak} --pheno ${phenotype} --grm \${grm_prefix} --reml \${phenotype_name}_LDAK-Thin --constrain YES
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
    // Extract "way" from VCF filename (all 5 types: SNP, INDEL, SV, SNP_INDEL, SNP_INDEL_SV)
    // VCF naming: SNP.split.xxx.vcf, INDEL.split.xxx.vcf, SV.split.xxx.vcf, SNP_INDEL.split.xxx.vcf, SNP_INDEL_SV.split.xxx.vcf
    // Note: Original VCFs are used for GEMMA. Recoding for LDAK happens internally in this workflow.
    gwas_with_way = gwas_input.map { vcf, pheno, list ->
        def vcf_name = vcf.name
        def way = ""
        if (vcf_name.startsWith("SNP_INDEL_SV.split.")) {
            way = "SNP_INDEL_SV_split"
        } else if (vcf_name.startsWith("SNP_INDEL.split.")) {
            way = "SNP_INDEL_split"
        } else if (vcf_name.startsWith("INDEL.split.")) {
            way = "INDEL_split"
        } else if (vcf_name.startsWith("SV.split.")) {
            way = "SV_split"
        } else if (vcf_name.startsWith("SNP.split.")) {
            way = "SNP_split"
        } else {
            way = "unknown_split"
        }
        tuple(vcf, pheno, list, way)
    }
    
    // Preprocess phenotype files
    preprocessed = PREPROCESS_PHENOTYPE_SPLIT(gwas_with_way)
    
    // Convert VCF to PLINK binary format
    plink_files = PLINK_VCF_CONVERSION(preprocessed)
    
    // Calculate kinship matrix with GEMMA
    kinship = GEMMA_KINSHIP(plink_files.plink_files)
    
    // Run GEMMA LMM association
    association = GEMMA_LMM_ASSOCIATION(kinship.kinship)
    
    // Perform LD clumping with PLINK
    clumped = PLINK_CLUMPING(association.association)
    
    // Estimate heritability with GCTA
    // heritability = GCTA_HERITABILITY(association.association)
    
    // Prepare LDAK GRM for each unique VCF+way combination
    unique_vcfs = gwas_with_way.map { vcf, pheno, list, way -> tuple(vcf, way) }.unique()
    
    // Recode VCFs for LDAK (only INDEL, SV, and merged VCFs - SNP doesn't need recoding)
    // Separate SNP VCFs from non-SNP VCFs
    snp_vcfs_for_ldak = unique_vcfs.filter { vcf, way -> way == "SNP_split" }
    non_snp_vcfs_for_ldak = unique_vcfs.filter { vcf, way -> way != "SNP_split" }
    
    // Recode non-SNP VCFs (INDEL, SV, and merged types)
    recoded_vcfs_for_ldak = RECODE_VCF_FOR_LDAK(non_snp_vcfs_for_ldak)
    
    // Combine recoded VCFs with SNP VCFs for LDAK GRM preparation
    all_vcfs_for_ldak = snp_vcfs_for_ldak.concat(recoded_vcfs_for_ldak.recoded_vcf)
    
    ldak_grm = LDAK_PREPARE_GRM(all_vcfs_for_ldak)
    
    // Combine GRM with phenotype information for LDAK heritability calculation
    // Reorganize both channels to put way at index 0 for reliable matching
    ldak_input = preprocessed
        .map { vcf, pheno, list, way -> tuple(way, vcf, pheno, list) }
        .combine(ldak_grm.grm_files, by: 0) // Match by way (index 0 in both)
        .map { way, vcf, pheno, list, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
            tuple(vcf, pheno, list, way, way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust) 
        }
    ldak_heritability = LDAK_CALCULATE_HERITABILITY(ldak_input)
    
    emit:
    assoc_files = association.association
    clumped_files = clumped.clumped
    // hsq_files = heritability.heritability
    ldak_reml_files = ldak_heritability.heritability
}

workflow GWAS_ANALYSIS_UNSPLIT {
    take:
    gwas_input  // tuple(vcf, phenotype, population_list)
    
    main:
    // Extract "way" from VCF filename (all 5 types: SNP, INDEL, SV, SNP_INDEL, SNP_INDEL_SV)
    // VCF naming: SNP.unsplit.xxx.vcf, INDEL.unsplit.xxx.vcf, SV.unsplit.xxx.vcf, SNP_INDEL.unsplit.xxx.vcf, SNP_INDEL_SV.unsplit.xxx.vcf
    // Note: Original VCFs are used for GEMMA. Recoding for LDAK happens internally in this workflow.
    gwas_with_way = gwas_input.map { vcf, pheno, list ->
        def vcf_name = vcf.name
        def way = ""
        if (vcf_name.startsWith("SNP_INDEL_SV.unsplit.")) {
            way = "SNP_INDEL_SV_unsplit"
        } else if (vcf_name.startsWith("SNP_INDEL.unsplit.")) {
            way = "SNP_INDEL_unsplit"
        } else if (vcf_name.startsWith("INDEL.unsplit.")) {
            way = "INDEL_unsplit"
        } else if (vcf_name.startsWith("SV.unsplit.")) {
            way = "SV_unsplit"
        } else if (vcf_name.startsWith("SNP.unsplit.")) {
            way = "SNP_unsplit"
        } else {
            way = "unknown_unsplit"
        }
        tuple(vcf, pheno, list, way)
    }
    
    // Preprocess phenotype files
    preprocessed = PREPROCESS_PHENOTYPE_UNSPLIT(gwas_with_way)
    
    // Convert VCF to PLINK binary format
    plink_files = PLINK_VCF_CONVERSION(preprocessed)
    
    // Calculate kinship matrix with GEMMA
    kinship = GEMMA_KINSHIP(plink_files.plink_files)
    
    // Run GEMMA LMM association
    association = GEMMA_LMM_ASSOCIATION(kinship.kinship)
    
    // Perform LD clumping with PLINK
    clumped = PLINK_CLUMPING(association.association)
    
    // Estimate heritability with GCTA
    // heritability = GCTA_HERITABILITY(association.association)
    
    // Prepare LDAK GRM for each unique VCF+way combination
    unique_vcfs = gwas_with_way.map { vcf, pheno, list, way -> tuple(vcf, way) }.unique()
    
    // Recode VCFs for LDAK (only INDEL, SV, and merged VCFs - SNP doesn't need recoding)
    // Separate SNP VCFs from non-SNP VCFs
    snp_vcfs_for_ldak = unique_vcfs.filter { vcf, way -> way == "SNP_unsplit" }
    non_snp_vcfs_for_ldak = unique_vcfs.filter { vcf, way -> way != "SNP_unsplit" }
    
    // Recode non-SNP VCFs (INDEL, SV, and merged types)
    recoded_vcfs_for_ldak = RECODE_VCF_FOR_LDAK(non_snp_vcfs_for_ldak)
    
    // Combine recoded VCFs with SNP VCFs for LDAK GRM preparation
    all_vcfs_for_ldak = snp_vcfs_for_ldak.concat(recoded_vcfs_for_ldak.recoded_vcf)
    
    ldak_grm = LDAK_PREPARE_GRM(all_vcfs_for_ldak)
    
    // Combine GRM with phenotype information for LDAK heritability calculation
    // Reorganize both channels to put way at index 0 for reliable matching
    ldak_input = preprocessed
        .map { vcf, pheno, list, way -> tuple(way, vcf, pheno, list) }
        .combine(ldak_grm.grm_files, by: 0) // Match by way (index 0 in both)
        .map { way, vcf, pheno, list, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
            tuple(vcf, pheno, list, way, way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust) 
        }
    ldak_heritability = LDAK_CALCULATE_HERITABILITY(ldak_input)
    
    emit:
    assoc_files = association.association
    clumped_files = clumped.clumped
    // hsq_files = heritability.heritability
    ldak_reml_files = ldak_heritability.heritability
}

