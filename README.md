# GraphPanGWAS - Nextflow GWAS Pipeline

A comprehensive Nextflow pipeline for genome-wide association studies (GWAS) of pangenome variants including SNPs, INDELs, and structural variants (SVs).

## Features

-  Automatic processing of both SPLIT and UNSPLIT variant analyses
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

## Some tips
### Ensure Your PanGenie SV VCF is Clean
Before running scripts, make sure your SV VCF file from PanGenie is clean and formatted properly.
For example:
```txt
P8.Nipponbare.TEJ.Chr1  1042    .       CC      C       .       PASS    UK=301;MA=14;ID=P8.Nipponbare.TEJ.Chr1-1042-DEL->4671>4673-1    GT      0|0     0|0
```
In particular, the last value in the INFO field should represent the SV length, as recommended in the official PanGenie guidelines: 
👉 PanGenie Wiki: [Genotyping variation nested inside of bubbles](https://github.com/eblerjana/pangenie/wiki/A:--Genotyping-variation-nested-inside-of-bubbles)
### Note on bcftools view
When using `bcftools view`, it may automatically add AC and AN fields to the INFO column, resulting in something like:
```txt
P8.Nipponbare.TEJ.Chr1  1042    .       CC      C       .       PASS    UK=301;MA=14;ID=P8.Nipponbare.TEJ.Chr1-1042-DEL->4671>4673-1;AC=20;AN=1396    GT
```
This modification can cause slient errors in scripts such as scripts/sv_filter.sh, which rely on the last number in the INFO field to represent SV length.
To avoid this issue, remove the automatically added AC and AN fields with:
```sh
bcftools annotate -x INFO/AC,INFO/AN input.vcf -o output.vcf
```

