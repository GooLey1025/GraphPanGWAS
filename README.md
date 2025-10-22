# GraphPanGWAS - Nextflow GWAS Pipeline

A comprehensive Nextflow pipeline for genome-wide association studies (GWAS) of pangenome variants including SNPs, INDELs, and structural variants (SVs).

## Features

- ✅Automatic processing of both SPLIT and UNSPLIT variant analyses
-  Parallel VCF preprocessing (SNP, INDEL, SV)
-  Automated GWAS with GEMMA (kinship-corrected LMM)
-  Lead marker identification with PLINK clumping
-  SNP-based heritability estimation with GCTA
-  Support for local and SLURM cluster execution
-  Automatic resume capability
-  Comprehensive execution reports

## Quick Start

```bash
# Run the complete pipeline
nextflow run main.nf -profile local -params-file 532rice.yaml -resume
```
