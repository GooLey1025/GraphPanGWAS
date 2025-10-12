#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    GraphPanGWAS - Nextflow GWAS Pipeline
========================================================================================
    Author: Converted from bash scripts
    Description: GWAS analysis pipeline for SNP, INDEL, and SV variants
========================================================================================
*/

// Print pipeline header
log.info """\
    ========================================
    G R A P H P A N G W A S   P I P E L I N E
    ========================================
    Population List    : ${params.population_list}
    SNP VCF            : ${params.snp_vcf}
    INDEL VCF          : ${params.indel_vcf}
    SV VCF             : ${params.sv_vcf}
    Phenotypes Dir     : ${params.phenotypes_dir}
    Output Prefix      : ${params.output_prefix}
    Base Analysis Type : ${params.analysis_type}
    Max Parallel GWAS  : ${params.max_parallel_gwas}
    
    Mode: BOTH split and unsplit analyses will be run
    ========================================
    """
    .stripIndent()

// Import modules
include { VCF_PREPROCESSING } from './modules/vcf_preprocessing'
include { GWAS_ANALYSIS } from './modules/gwas'
include { EXTRACT_HERITABILITY } from './modules/postprocessing'
include { EXTRACT_LEAD_MARKERS } from './modules/postprocessing'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Create channels for input files
    snp_vcf_ch = params.snp_vcf ? Channel.fromPath(params.snp_vcf, checkIfExists: true) : Channel.empty()
    indel_vcf_ch = params.indel_vcf ? Channel.fromPath(params.indel_vcf, checkIfExists: true) : Channel.empty()
    sv_vcf_ch = params.sv_vcf ? Channel.fromPath(params.sv_vcf, checkIfExists: true) : Channel.empty()
    population_list_ch = Channel.fromPath(params.population_list, checkIfExists: true)
    
    // VCF Preprocessing based on analysis type
    VCF_PREPROCESSING(
        snp_vcf_ch,
        indel_vcf_ch,
        sv_vcf_ch,
        params.analysis_type
    )
    
    // Validate and collect phenotype files
    phenotype_files = Channel
        .fromPath("${params.phenotypes_dir}/*.tsv", checkIfExists: true)
        .map { file ->
            def lines = file.readLines()
            if (lines.size() < 2) {
                log.error "ERROR: ${file.name} has less than 2 lines"
                System.exit(1)
            }
            def header = lines[0].split('\t')
            if (header.size() != 2) {
                log.error "ERROR: ${file.name} does not have exactly 2 columns (has ${header.size()} columns)"
                log.error "       Expected format: ID<tab>phenotype_value"
                System.exit(1)
            }
            if (header[0] != 'ID') {
                log.error "ERROR: ${file.name} first column is '${header[0]}', expected 'ID'"
                log.error "       Header line should be: ID<tab>phenotype_name"
                System.exit(1)
            }
            return file
        }
    
    // ===== Run GWAS on SPLIT version =====
    gwas_input_split = VCF_PREPROCESSING.out.split_vcf
        .combine(phenotype_files)
        .combine(population_list_ch)
    
    GWAS_ANALYSIS_SPLIT: GWAS_ANALYSIS(
        gwas_input_split.map { vcf, pheno, list -> tuple(vcf, pheno, list) }
    )
    
    // Post-processing for SPLIT
    EXTRACT_HERITABILITY(
        GWAS_ANALYSIS_SPLIT.out.hsq_files.collect(),
        params.output_prefix,
        "${params.analysis_type}_split"
    )
    
    EXTRACT_LEAD_MARKERS(
        GWAS_ANALYSIS_SPLIT.out.clumped_files.collect(),
        params.output_prefix,
        "${params.analysis_type}_split"
    )
    
    // ===== Run GWAS on UNSPLIT version =====
    gwas_input_unsplit = VCF_PREPROCESSING.out.unsplit_vcf
        .combine(phenotype_files)
        .combine(population_list_ch)
    
    GWAS_ANALYSIS_UNSPLIT: GWAS_ANALYSIS(
        gwas_input_unsplit.map { vcf, pheno, list -> tuple(vcf, pheno, list) }
    )
    
    // Post-processing for UNSPLIT
    EXTRACT_HERITABILITY(
        GWAS_ANALYSIS_UNSPLIT.out.hsq_files.collect(),
        params.output_prefix,
        "${params.analysis_type}_unsplit"
    )
    
    EXTRACT_LEAD_MARKERS(
        GWAS_ANALYSIS_UNSPLIT.out.clumped_files.collect(),
        params.output_prefix,
        "${params.analysis_type}_unsplit"
    )
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        ========================================
        Pipeline completed!
        Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration  : ${workflow.duration}
        
        Results:
          SPLIT   : ${params.output_prefix}/${params.analysis_type}_split/
          UNSPLIT : ${params.output_prefix}/${params.analysis_type}_unsplit/
        ========================================
        """
        .stripIndent()
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}

