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
    Max Parallel GWAS  : ${params.max_parallel_gwas}
    
    Mode: BOTH split and unsplit analyses will be run
    ========================================
    """
    .stripIndent()

// Import modules
include { VCF_PREPROCESSING } from './modules/vcf_preprocessing'
include { WHOLE_GENOME_MARKERS } from './modules/wholegenome_markers'
include { GWAS_ANALYSIS_SPLIT } from './modules/gwas'
include { GWAS_ANALYSIS_UNSPLIT } from './modules/gwas'
// ===== GCTA HERITABILITY EXTRACTION MODULES - COMMENTED OUT =====
// include { EXTRACT_HERITABILITY as EXTRACT_HERITABILITY_SPLIT } from './modules/postprocessing'
// include { EXTRACT_HERITABILITY as EXTRACT_HERITABILITY_UNSPLIT } from './modules/postprocessing'
include { EXTRACT_LDAK_HERITABILITY as EXTRACT_LDAK_HERITABILITY_SPLIT } from './modules/postprocessing'
include { GENERATE_HERITABILITY_TABLE } from './modules/postprocessing'
include { EXTRACT_LEAD_MARKERS as EXTRACT_LEAD_MARKERS_SPLIT } from './modules/postprocessing'
include { EXTRACT_LEAD_MARKERS as EXTRACT_LEAD_MARKERS_UNSPLIT } from './modules/postprocessing'
include { EXTRACT_LD_MARKER_VCF } from './modules/postprocessing'

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
    ref_fa_ch = Channel.fromPath(params.ref_fa, checkIfExists: true)
    // VCF Preprocessing (generates both split and unsplit for all available types)
    VCF_PREPROCESSING(
        snp_vcf_ch,
        indel_vcf_ch,
        sv_vcf_ch,
        ref_fa_ch
    )
    
    // Generate whole genome LD-pruned markers for SNP split VCF only
    // Use the SNP split VCF from VCF preprocessing output
    snp_split_vcf = VCF_PREPROCESSING.out.split_vcf
        .filter { vcf -> vcf.name.startsWith("SNP.split.") }
    
    // Run whole genome marker generation
    whole_genome_markers = WHOLE_GENOME_MARKERS(snp_split_vcf, population_list_ch)
    
    // Validate and collect phenotype files
    phenotype_files = Channel
        .fromPath("${params.phenotypes_dir}/*.tsv", checkIfExists: true)
        .map { file ->
            def lines = file.readLines()
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
    
    gwas_results_split = GWAS_ANALYSIS_SPLIT(
        gwas_input_split.map { vcf, pheno, list -> tuple(vcf, pheno, list) }
    )
    
    // Post-processing for SPLIT - group by way
    // hsq_by_way_split = gwas_results_split.hsq_files
    //     .map { pheno_name, hsq, way -> tuple(way, hsq) }
    //     .groupTuple()
    
    clumped_by_way_split = gwas_results_split.clumped_files
        .map { pheno_name, clumped, way -> tuple(way, clumped) }
        .groupTuple()
    
    EXTRACT_LEAD_MARKERS_SPLIT(
        clumped_by_way_split,
        params.output_prefix
    )
    
    // ===== Extract LD marker-filtered VCFs for split version =====
    // Only for SNP_split, SNP_INDEL_split, and SNP_INDEL_SV_split
    ld_marker_files = EXTRACT_LEAD_MARKERS_SPLIT.out.lead_markers_summary
    
    // Get the original split VCFs for the three ways
    split_vcfs_for_ld = VCF_PREPROCESSING.out.split_vcf
        .map { vcf ->
            def vcf_name = vcf.name
            def way = ""
            if (vcf_name.startsWith("SNP_INDEL_SV.split.")) {
                way = "SNP_INDEL_SV_split"
            } else if (vcf_name.startsWith("SNP_INDEL.split.")) {
                way = "SNP_INDEL_split"
            } else if (vcf_name.startsWith("SNP.split.")) {
                way = "SNP_split"
            }
            if (way != "") {
                tuple(way, vcf)
            }
        }
        .filter { it != null }
    
    // Combine LD marker files with their corresponding VCFs and population list
    ld_extraction_input = ld_marker_files.combine(split_vcfs_for_ld, by: 0).combine(population_list_ch)
    
    EXTRACT_LD_MARKER_VCF(ld_extraction_input)
    
    // ===== Run GWAS on UNSPLIT version =====
    gwas_input_unsplit = VCF_PREPROCESSING.out.unsplit_vcf
        .combine(phenotype_files)
        .combine(population_list_ch)
    
    gwas_results_unsplit = GWAS_ANALYSIS_UNSPLIT(
        gwas_input_unsplit.map { vcf, pheno, list -> tuple(vcf, pheno, list) }
    )
    
    // Post-processing for UNSPLIT - group by way
    // hsq_by_way_unsplit = gwas_results_unsplit.hsq_files
    //     .map { pheno_name, hsq, way -> tuple(way, hsq) }
    //     .groupTuple()
    
    clumped_by_way_unsplit = gwas_results_unsplit.clumped_files
        .map { pheno_name, clumped, way -> tuple(way, clumped) }
        .groupTuple()
    
    EXTRACT_LEAD_MARKERS_UNSPLIT(
        clumped_by_way_unsplit,
        params.output_prefix
    )
    
    // ===== Combine all LDAK heritability results (all 5 types: SNP, INDEL, SV, SNP_INDEL, SNP_INDEL_SV) =====
    // Merge split and unsplit LDAK results together
    all_ldak_reml = gwas_results_split.ldak_reml_files
        .concat(gwas_results_unsplit.ldak_reml_files)
        .map { pheno_name, reml, way -> tuple(way, reml) }
        .groupTuple()
    
    // Extract LDAK heritability for all types together (both split and unsplit)
    EXTRACT_LDAK_HERITABILITY_SPLIT(
        all_ldak_reml,
        params.output_prefix
    )
    
    // Collect all LDAK TSV files and generate Excel summary table
    all_ldak_tsv_files = EXTRACT_LDAK_HERITABILITY_SPLIT.out.ldak_heritability_summary
        .collect()
    
    GENERATE_HERITABILITY_TABLE(
        all_ldak_tsv_files,
        params.output_prefix,
        params.heritability_dir
    )
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

// workflow.onComplete {
//     log.info """\
//         ========================================
//         Pipeline completed!
//         Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
//         Duration  : ${workflow.duration}
        
//         Results:
//           SPLIT   : ${params.output_prefix}/${analysis_type_name}_split/
//           UNSPLIT : ${params.output_prefix}/${analysis_type_name}_unsplit/
//         ========================================
//         """
//         .stripIndent()
// }

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}

