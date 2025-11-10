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
        --double-id \\
        --memory ${params.plink_memory} \\
        --pheno ${phenotype} \\
        --mpheno 1 \\
        --make-bed \\
        --out ${phenotype.baseName} \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --maf ${params.maf_threshold}
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

process MANHATTAN_PLOT {
    tag "${phenotype_name}_${way}"
    label 'process_low'
    publishDir "${params.output_prefix}/${way}/plots/", 
               mode: 'copy', 
               pattern: "*.manhattan.png"
    
    input:
    tuple val(phenotype_name), path(bed), path(bim), path(fam), path(assoc), path(phenotype), val(way)
    
    output:
    path("${phenotype_name}.manhattan.png"), emit: plot
    
    script:
    """
    ${params.python3} ${projectDir}/scripts/manhattan_plot.py \\
        ${assoc} \\
        ${phenotype_name}.manhattan.png \\
        processed_${phenotype_name}_${way}
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
    plink --vcf ${vcf} --double-id --memory 10000 --make-bed --out ${vcf.baseName}_WG --allow-extra-chr --allow-no-sex --maf ${params.maf_threshold}

    # Calculate GRM with user-specified power parameter
    ${params.ldak} --calc-kins-direct ${vcf.baseName}_WG_LDAK-Thin --bfile ${vcf.baseName}_WG --power ${params.ldak_power}
    """
}

process LDAK_BATCH_HERITABILITY_SPLIT {
    tag "LDAK_BATCH_split"
    label 'process_high'
    publishDir "${params.output_prefix}", mode: 'copy', pattern: "*_split/results/**/*.reml"
    
    input:
    tuple path(snp_grm_bin, stageAs: 'snp.grm.bin'), 
          path(snp_grm_id, stageAs: 'snp.grm.id'), 
          path(snp_grm_details, stageAs: 'snp.grm.details'), 
          path(snp_grm_adjust, stageAs: 'snp.grm.adjust'),
          path(indel_grm_bin, stageAs: 'indel.grm.bin'), 
          path(indel_grm_id, stageAs: 'indel.grm.id'), 
          path(indel_grm_details, stageAs: 'indel.grm.details'), 
          path(indel_grm_adjust, stageAs: 'indel.grm.adjust'),
          path(sv_grm_bin, stageAs: 'sv.grm.bin'), 
          path(sv_grm_id, stageAs: 'sv.grm.id'), 
          path(sv_grm_details, stageAs: 'sv.grm.details'), 
          path(sv_grm_adjust, stageAs: 'sv.grm.adjust'),
          val(phenotypes)
    
    output:
    path("*_split/results/**/*.reml"), emit: reml_files
    
    script:
    // phenotypes is a string with file paths separated by newlines
    // Parse it back to a List
    def pheno_paths_str = phenotypes.toString().trim()
    def pheno_paths = pheno_paths_str ? pheno_paths_str.split('\n') as List : []
    
    def pheno_files = pheno_paths.collect { path_str ->
        def f = file(path_str.trim())
        f.toAbsolutePath().toString()
    }
    def pheno_files_str = pheno_files.join(' ')
    def pheno_count = pheno_files.size()
    """
    # Get GRM prefixes (now using staged names)
    snp_grm="snp"
    indel_grm="indel"
    sv_grm="sv"
    
    # Create GRM list files for merged types
    printf "\${snp_grm}\\n\${indel_grm}\\n" > SNP_INDEL.grm.list
    printf "\${snp_grm}\\n\${indel_grm}\\n\${sv_grm}\\n" > SNP_INDEL_SV.grm.list
    
    # Create phenotypes directory and copy all phenotype files
    mkdir -p phenotypes
    echo "=========================================="
    echo "Total phenotype files to copy: ${pheno_count}"
    echo "Phenotype file paths:"
    echo "${pheno_files_str}"
    echo "=========================================="
    
    # Copy each phenotype file with detailed logging
    copy_count=0
    for pheno_file in ${pheno_files_str}; do
        echo "Processing: \${pheno_file}"
        if [ -f "\${pheno_file}" ]; then
            filename=\$(basename "\${pheno_file}")
            cp "\${pheno_file}" "phenotypes/\${filename}"
            if [ \$? -eq 0 ]; then
                echo "  SUCCESS: Copied \${pheno_file} -> phenotypes/\${filename}"
                copy_count=\$((copy_count + 1))
            else
                echo "  ERROR: Failed to copy \${pheno_file}"
                exit 1
            fi
        else
            echo "  ERROR: File not found: \${pheno_file}"
            echo "  Current directory: \$(pwd)"
            echo "  File exists check: [ -f \"\${pheno_file}\" ]"
            exit 1
        fi
    done
    
    echo "=========================================="
    echo "Successfully copied \${copy_count} phenotype files"
    echo "Verifying copied files in phenotypes directory:"
    ls -lh phenotypes/ | head -30
    echo "Total files in phenotypes directory:"
    ls phenotypes/* 2>/dev/null | wc -l || echo "0"
    echo "=========================================="
    
    # Create output directories
    mkdir -p SNP_split/results INDEL_split/results SV_split/results
    mkdir -p SNP_INDEL_split/results SNP_INDEL_SV_split/results
    
    # Process all phenotypes in parallel
    export snp_grm indel_grm sv_grm
    export LDAK="${params.ldak}"
    
    # Check if phenotypes directory has files before processing
    if [ -n "\$(ls -A phenotypes/* 2>/dev/null)" ]; then
        ls phenotypes/* | parallel -j ${task.cpus} '
        pheno_file={}
        pheno_name=\$(basename {} | sed "s/\\.[^.]*\$//")
        
        # Create directories for this phenotype
        mkdir -p SNP_split/results/\${pheno_name}
        mkdir -p INDEL_split/results/\${pheno_name}
        mkdir -p SV_split/results/\${pheno_name}
        mkdir -p SNP_INDEL_split/results/\${pheno_name}
        mkdir -p SNP_INDEL_SV_split/results/\${pheno_name}
        
        # SNP heritability
        \$LDAK --reml SNP_split/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --grm \$snp_grm --pheno \$pheno_file --constrain YES
        
        # INDEL heritability
        \$LDAK --reml INDEL_split/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --grm \$indel_grm --pheno \$pheno_file --constrain YES
        
        # SV heritability
        \$LDAK --reml SV_split/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --grm \$sv_grm --pheno \$pheno_file --constrain YES
        
        # SNP_INDEL merged heritability
        \$LDAK --reml SNP_INDEL_split/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --mgrm SNP_INDEL.grm.list --pheno \$pheno_file --constrain YES
        
        # SNP_INDEL_SV merged heritability
        \$LDAK --reml SNP_INDEL_SV_split/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --mgrm SNP_INDEL_SV.grm.list --pheno \$pheno_file --constrain YES
        '
    else
        echo "ERROR: No phenotype files found in phenotypes directory!"
        exit 1
    fi
    """
}

process LDAK_BATCH_HERITABILITY_UNSPLIT {
    tag "LDAK_BATCH_unsplit"
    label 'process_high'
    publishDir "${params.output_prefix}", mode: 'copy', pattern: "*_unsplit/results/**/*.reml"
    
    input:
    tuple path(snp_grm_bin, stageAs: 'snp.grm.bin'), 
          path(snp_grm_id, stageAs: 'snp.grm.id'), 
          path(snp_grm_details, stageAs: 'snp.grm.details'), 
          path(snp_grm_adjust, stageAs: 'snp.grm.adjust'),
          path(indel_grm_bin, stageAs: 'indel.grm.bin'), 
          path(indel_grm_id, stageAs: 'indel.grm.id'), 
          path(indel_grm_details, stageAs: 'indel.grm.details'), 
          path(indel_grm_adjust, stageAs: 'indel.grm.adjust'),
          path(sv_grm_bin, stageAs: 'sv.grm.bin'), 
          path(sv_grm_id, stageAs: 'sv.grm.id'), 
          path(sv_grm_details, stageAs: 'sv.grm.details'), 
          path(sv_grm_adjust, stageAs: 'sv.grm.adjust'),
          val(phenotypes)
    
    output:
    path("*_unsplit/results/**/*.reml"), emit: reml_files
    
    script:
    // phenotypes is a string with file paths separated by newlines
    // Parse it back to a List
    def pheno_paths_str = phenotypes.toString().trim()
    def pheno_paths = pheno_paths_str ? pheno_paths_str.split('\n') as List : []
    
    def pheno_files = pheno_paths.collect { path_str ->
        def f = file(path_str.trim())
        f.toAbsolutePath().toString()
    }
    def pheno_files_str = pheno_files.join(' ')
    def pheno_count = pheno_files.size()
    """
    # Get GRM prefixes (now using staged names)
    snp_grm="snp"
    indel_grm="indel"
    sv_grm="sv"
    
    # Create GRM list files for merged types
    printf "\${snp_grm}\\n\${indel_grm}\\n" > SNP_INDEL.grm.list
    printf "\${snp_grm}\\n\${indel_grm}\\n\${sv_grm}\\n" > SNP_INDEL_SV.grm.list
    
    # Create phenotypes directory and copy all phenotype files
    mkdir -p phenotypes
    echo "=========================================="
    echo "Total phenotype files to copy: ${pheno_count}"
    echo "Phenotype file paths:"
    echo "${pheno_files_str}"
    echo "=========================================="
    
    # Copy each phenotype file with detailed logging
    copy_count=0
    for pheno_file in ${pheno_files_str}; do
        echo "Processing: \${pheno_file}"
        if [ -f "\${pheno_file}" ]; then
            filename=\$(basename "\${pheno_file}")
            cp "\${pheno_file}" "phenotypes/\${filename}"
            if [ \$? -eq 0 ]; then
                echo "  SUCCESS: Copied \${pheno_file} -> phenotypes/\${filename}"
                copy_count=\$((copy_count + 1))
            else
                echo "  ERROR: Failed to copy \${pheno_file}"
                exit 1
            fi
        else
            echo "  ERROR: File not found: \${pheno_file}"
            echo "  Current directory: \$(pwd)"
            echo "  File exists check: [ -f \"\${pheno_file}\" ]"
            exit 1
        fi
    done
    
    echo "=========================================="
    echo "Successfully copied \${copy_count} phenotype files"
    echo "Verifying copied files in phenotypes directory:"
    ls -lh phenotypes/ | head -30
    echo "Total files in phenotypes directory:"
    ls phenotypes/* 2>/dev/null | wc -l || echo "0"
    echo "=========================================="
    
    # Create output directories
    mkdir -p SNP_unsplit/results INDEL_unsplit/results SV_unsplit/results
    mkdir -p SNP_INDEL_unsplit/results SNP_INDEL_SV_unsplit/results
    
    # Process all phenotypes in parallel
    export snp_grm indel_grm sv_grm
    export LDAK="${params.ldak}"
    
    # Check if phenotypes directory has files before processing
    if [ -n "\$(ls -A phenotypes/* 2>/dev/null)" ]; then
        ls phenotypes/* | parallel -j ${task.cpus} '
        pheno_file={}
        pheno_name=\$(basename {} | sed "s/\\.[^.]*\$//")
        
        # Create directories for this phenotype
        mkdir -p SNP_unsplit/results/\${pheno_name}
        mkdir -p INDEL_unsplit/results/\${pheno_name}
        mkdir -p SV_unsplit/results/\${pheno_name}
        mkdir -p SNP_INDEL_unsplit/results/\${pheno_name}
        mkdir -p SNP_INDEL_SV_unsplit/results/\${pheno_name}
        
        # SNP heritability
        \$LDAK --reml SNP_unsplit/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --grm \$snp_grm --pheno \$pheno_file --constrain YES
        
        # INDEL heritability
        \$LDAK --reml INDEL_unsplit/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --grm \$indel_grm --pheno \$pheno_file --constrain YES
        
        # SV heritability
        \$LDAK --reml SV_unsplit/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --grm \$sv_grm --pheno \$pheno_file --constrain YES
        
        # SNP_INDEL merged heritability
        \$LDAK --reml SNP_INDEL_unsplit/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --mgrm SNP_INDEL.grm.list --pheno \$pheno_file --constrain YES
        
        # SNP_INDEL_SV merged heritability
        \$LDAK --reml SNP_INDEL_SV_unsplit/results/\${pheno_name}/\${pheno_name}_LDAK-Thin \\
            --mgrm SNP_INDEL_SV.grm.list --pheno \$pheno_file --constrain YES
        '
    else
        echo "ERROR: No phenotype files found in phenotypes directory!"
        exit 1
    fi
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
    
    // plot the manhattan plot of the association results
    MANHATTAN_PLOT(association.association)
    
    // Perform LD clumping with PLINK
    clumped = PLINK_CLUMPING(association.association)
    
    // Estimate heritability with GCTA
    // heritability = GCTA_HERITABILITY(association.association)
    
    // Prepare LDAK GRM for individual variant types (SNP, INDEL, SV)
    individual_vcfs = gwas_with_way
        .map { vcf, pheno, list, way -> tuple(vcf, way) }
        .unique()
        .filter { vcf, way -> 
            way in ["SNP_split", "INDEL_split", "SV_split"]
        }
    
    // Recode all VCFs for LDAK
    recoded_vcfs_for_ldak = RECODE_VCF_FOR_LDAK(individual_vcfs)
    
    // Generate GRMs for individual variant types
    ldak_grm = LDAK_PREPARE_GRM(recoded_vcfs_for_ldak.recoded_vcf)
    
    // Collect all 3 GRMs and all phenotypes for batch processing
    snp_grm = ldak_grm.grm_files.filter { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
        way == "SNP_split" 
    }.map { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(grm_bin, grm_id, grm_details, grm_adjust)
    }
    
    indel_grm = ldak_grm.grm_files.filter { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
        way == "INDEL_split" 
    }.map { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(grm_bin, grm_id, grm_details, grm_adjust)
    }
    
    sv_grm = ldak_grm.grm_files.filter { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
        way == "SV_split" 
    }.map { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(grm_bin, grm_id, grm_details, grm_adjust)
    }
    
    // Collect all phenotypes (unique by file name to avoid duplicates across different VCFs)
    all_phenotypes = preprocessed
        .map { vcf, pheno, list, way -> tuple(pheno.name, pheno) }
        .unique { it[0] }
        .map { name, pheno -> pheno }
        .collect()
    
    // Combine all GRMs and phenotypes for batch processing
    // Problem: .combine() expands Lists, so we need to avoid using combine with the List
    // Solution: Convert the List to a string (file paths separated by newlines), then parse in process
    grm_tuple = snp_grm.first()
        .combine(indel_grm.first())
        .combine(sv_grm.first())
    
    // Convert phenotype list to a string of file paths (one per line)
    // This prevents .combine() from expanding the List
    pheno_list_str = all_phenotypes
        .map { pheno_list ->
            def paths = pheno_list.collect { it.toString() }
            paths.join('\n')
        }
    
    // Combine GRM tuple with phenotype list string
    batch_input = grm_tuple
        .combine(pheno_list_str)
        .map { snp_bin, snp_id, snp_details, snp_adjust, indel_bin, indel_id, indel_details, indel_adjust, sv_bin, sv_id, sv_details, sv_adjust, pheno_str ->
            // Parse the string back to a List in the process
            tuple(snp_bin, snp_id, snp_details, snp_adjust,
                  indel_bin, indel_id, indel_details, indel_adjust,
                  sv_bin, sv_id, sv_details, sv_adjust,
                  pheno_str)
        }
    
    // Run batch heritability calculation
    ldak_heritability_raw = LDAK_BATCH_HERITABILITY_SPLIT(batch_input)
    
    // Transform output to match expected format: (pheno_name, reml, way)
    ldak_heritability = ldak_heritability_raw.reml_files
        .flatten()
        .map { reml_file ->
            // Extract way and phenotype from path like "SNP_split/results/phenotype_name/phenotype_name_LDAK-Thin.reml"
            // Use file.name to get just the filename, and navigate up the directory structure
            def file_path = reml_file.toString()
            def path_parts = file_path.tokenize('/')
            // Find the index of "results" and extract way from the part before it
            def results_idx = path_parts.findIndexOf { it == 'results' }
            def way = results_idx > 0 ? path_parts[results_idx - 1] : path_parts[-4]  // e.g., "SNP_split"
            def pheno_name = path_parts[results_idx + 1]  // e.g., "phenotype_name"
            tuple(pheno_name, reml_file, way)
        }
    
    emit:
    assoc_files = association.association
    clumped_files = clumped.clumped
    ldak_reml_files = ldak_heritability
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
    
    // plot the manhattan plot of the association results
    MANHATTAN_PLOT(association.association)
    
    // Perform LD clumping with PLINK
    clumped = PLINK_CLUMPING(association.association)
    
    // Estimate heritability with GCTA
    // heritability = GCTA_HERITABILITY(association.association)
    
    // Prepare LDAK GRM for individual variant types (SNP, INDEL, SV)
    individual_vcfs = gwas_with_way
        .map { vcf, pheno, list, way -> tuple(vcf, way) }
        .unique()
        .filter { vcf, way -> 
            way in ["SNP_unsplit", "INDEL_unsplit", "SV_unsplit"]
        }
    
    // Recode all VCFs for LDAK
    recoded_vcfs_for_ldak = RECODE_VCF_FOR_LDAK(individual_vcfs)
    
    // Generate GRMs for individual variant types
    ldak_grm = LDAK_PREPARE_GRM(recoded_vcfs_for_ldak.recoded_vcf)
    
    // Collect all 3 GRMs and all phenotypes for batch processing
    snp_grm = ldak_grm.grm_files.filter { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
        way == "SNP_unsplit" 
    }.map { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(grm_bin, grm_id, grm_details, grm_adjust)
    }
    
    indel_grm = ldak_grm.grm_files.filter { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
        way == "INDEL_unsplit" 
    }.map { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(grm_bin, grm_id, grm_details, grm_adjust)
    }
    
    sv_grm = ldak_grm.grm_files.filter { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust -> 
        way == "SV_unsplit" 
    }.map { way, bed, bim, fam, grm_bin, grm_id, grm_details, grm_adjust ->
        tuple(grm_bin, grm_id, grm_details, grm_adjust)
    }
    
    // Collect all phenotypes (unique by file name to avoid duplicates across different VCFs)
    all_phenotypes = preprocessed
        .map { vcf, pheno, list, way -> tuple(pheno.name, pheno) }
        .unique { it[0] }
        .map { name, pheno -> pheno }
        .collect()
    
    // Combine all GRMs and phenotypes for batch processing
    // Problem: .combine() expands Lists, so we need to avoid using combine with the List
    // Solution: Convert the List to a string (file paths separated by newlines), then parse in process
    grm_tuple = snp_grm.first()
        .combine(indel_grm.first())
        .combine(sv_grm.first())
    
    // Convert phenotype list to a string of file paths (one per line)
    // This prevents .combine() from expanding the List
    pheno_list_str = all_phenotypes
        .map { pheno_list ->
            def paths = pheno_list.collect { it.toString() }
            paths.join('\n')
        }
    
    // Combine GRM tuple with phenotype list string
    batch_input = grm_tuple
        .combine(pheno_list_str)
        .map { snp_bin, snp_id, snp_details, snp_adjust, indel_bin, indel_id, indel_details, indel_adjust, sv_bin, sv_id, sv_details, sv_adjust, pheno_str ->
            // Parse the string back to a List in the process
            tuple(snp_bin, snp_id, snp_details, snp_adjust,
                  indel_bin, indel_id, indel_details, indel_adjust,
                  sv_bin, sv_id, sv_details, sv_adjust,
                  pheno_str)
        }
    
    // Run batch heritability calculation
    ldak_heritability_raw = LDAK_BATCH_HERITABILITY_UNSPLIT(batch_input)
    
    // Transform output to match expected format: (pheno_name, reml, way)
    ldak_heritability = ldak_heritability_raw.reml_files
        .flatten()
        .map { reml_file ->
            // Extract way and phenotype from path like "SNP_unsplit/results/phenotype_name/phenotype_name_LDAK-Thin.reml"
            // Use file.name to get just the filename, and navigate up the directory structure
            def file_path = reml_file.toString()
            def path_parts = file_path.tokenize('/')
            // Find the index of "results" and extract way from the part before it
            def results_idx = path_parts.findIndexOf { it == 'results' }
            def way = results_idx > 0 ? path_parts[results_idx - 1] : path_parts[-4]  // e.g., "SNP_unsplit"
            def pheno_name = path_parts[results_idx + 1]  // e.g., "phenotype_name"
            tuple(pheno_name, reml_file, way)
        }
    
    emit:
    assoc_files = association.association
    clumped_files = clumped.clumped
    ldak_reml_files = ldak_heritability
}

