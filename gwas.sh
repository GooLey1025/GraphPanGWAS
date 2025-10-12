#!/bin/bash
set -euo pipefail

# Parse arguments
phenotype_csv=$1
population_list=$2
vcf=$3
P=$4
way=$5
phenotype=$(basename "$phenotype_csv")

# Tool paths
gemma=~/biosoftwares/GEMMA/gemma-0.98.5-linux-static-AMD64

# Logging functions - output to stdout/stderr (will be redirected by caller)
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1: $2"; }
log_info() { log "INFO" "$1"; }
log_error() { log "ERROR" "$1" >&2; }

log_info "Starting GWAS: phenotype=$phenotype, vcf=$vcf, way=$way"

# Initialize working directory
work_dir="$P/$way"
mkdir -p "$work_dir" || { log_error "Failed to create $work_dir"; exit 1; }
cd "$work_dir" || { log_error "Failed to cd to $work_dir"; exit 1; }

mkdir -p tmp phes_tmp gemma_tmp "tmp/$phenotype" "results/$phenotype" lead_markers gcta_tmp || {
  log_error "Failed to create working directories"; exit 1;
}
export TMPDIR=./tmp

# Phenotype preprocessing
log_info "Preprocessing phenotype"
python3 ../../scripts/process_phenotype.py "../../$population_list" "../../$phenotype_csv" "./phes_tmp/sorted_$phenotype" || {
  log_error "Phenotype processing failed"; exit 1;
}

awk '{printf "%s", $1; for(i=1; i<=NF; i++) printf "\t%s", $i; print "";}' \
  "./phes_tmp/sorted_$phenotype" > "./phes_tmp/tmp_$phenotype.pre_plink.csv" || {
  log_error "Failed to format phenotype"; exit 1;
}

awk 'BEGIN {FS="\t"; OFS="\t"} NR==1 {$1="FID"; $2="IID"} {print}' \
  "./phes_tmp/tmp_$phenotype.pre_plink.csv" > "./phes_tmp/$phenotype.pre_plink.csv" || {
  log_error "Failed to add FID/IID headers"; exit 1;
}

# PLINK: Convert VCF to binary format
log_info "Running PLINK VCF conversion"
plink --vcf "../../vcf_dir/$vcf" --memory 10000 --pheno "./phes_tmp/$phenotype.pre_plink.csv" \
  --mpheno 1 --make-bed --out "./phes_tmp/$phenotype" --allow-extra-chr --allow-no-sex || {
  log_error "PLINK conversion failed"; exit 1;
}

# GEMMA: Kinship matrix and LMM association
log_info "Running GEMMA kinship calculation"
$gemma -bfile "./phes_tmp/$phenotype" -gk 1 -outdir "./tmp/$phenotype/" -o kinship \
  -miss 0.2 -maf 0.01 -hwe 0 || { log_error "GEMMA kinship failed"; exit 1; }

log_info "Running GEMMA LMM association"
$gemma -bfile "./phes_tmp/$phenotype" -k "./tmp/$phenotype/kinship.cXX.txt" \
  -outdir "results/$phenotype" -o gemma_lmm -lmm 1 -miss 0.2 -maf 0.01 || {
  log_error "GEMMA LMM failed"; exit 1;
}

# PLINK: Lead marker clumping
log_info "Running PLINK clumping"
plink --bfile "./phes_tmp/$phenotype" --clump "results/$phenotype/gemma_lmm.assoc.txt" \
  --clump-p1 1e-6 --clump-p2 0.05 --clump-r2 0.1 --clump-kb 1000 --allow-no-sex \
  --out "lead_markers/$phenotype" --clump-snp-field rs --clump-field p_wald || {
  log_error "PLINK clumping failed"; exit 1;
}

# GCTA: Heritability estimation
log_info "Running GCTA heritability analysis"
cd gcta_tmp || { log_error "Failed to cd to gcta_tmp"; exit 1; }

gcta64 --bfile "../phes_tmp/$phenotype" --autosome --maf 0.01 --make-grm \
  --out "$phenotype" --thread-num 10 || { log_error "GCTA GRM failed"; exit 1; }

gcta64 --grm "$phenotype" --pheno "../phes_tmp/$phenotype.pre_plink.csv" --reml \
  --out "../results/$phenotype/$phenotype" --thread-num 10 || {
  log_error "GCTA REML failed"; exit 1;
}

log_info "GWAS completed successfully"

