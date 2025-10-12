# GraphPanGWAS - Nextflow Pipeline

A Nextflow workflow for genome-wide association studies (GWAS) of pangenome variants including SNPs, INDELs, and structural variants (SVs).

## Overview

This pipeline performs comprehensive GWAS analysis including:

1. **VCF Preprocessing**: Chromosome renaming, variant normalization, ID assignment, and merging
2. **GWAS Analysis**: Association testing using GEMMA's linear mixed models with kinship correction
3. **Lead Marker Identification**: LD-based clumping with PLINK
4. **Heritability Estimation**: Variance component estimation with GCTA
5. **Results Summarization**: Extraction and compilation of key results

## Requirements

### Software Dependencies

The pipeline requires the following tools to be installed and available in your PATH (or specify paths in configuration):

- **Nextflow** >= 22.10.0
- **PLINK** (1.9 or 2.0)
- **GEMMA** (>= 0.98)
- **GCTA** (>= 1.93)
- **bcftools** (>= 1.10)
- **bgzip** (tabix package)
- **Python 3** with pandas
- **awk**, **perl** (usually pre-installed on Linux)

### Installation

#### Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/  # or add to your PATH
```

#### Install Analysis Tools

```bash
# Example using conda
conda install -c bioconda plink bcftools tabix

# GEMMA - download from https://github.com/genetics-statistics/GEMMA/releases
wget https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
gunzip gemma-0.98.5-linux-static-AMD64.gz
chmod +x gemma-0.98.5-linux-static-AMD64
mv gemma-0.98.5-linux-static-AMD64 ~/bin/gemma

# GCTA - download from https://yanglab.westlake.edu.cn/software/gcta/
wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip
unzip gcta-1.94.1-linux-kernel-3-x86_64.zip
mv gcta-1.94.1-linux-kernel-3-x86_64/gcta64 ~/bin/
```

## Input Files

### Required Inputs

1. **Population List** (`--population_list`): Tab-delimited file with sample IDs
   ```
   Sample001
   Sample002
   Sample003
   ```

2. **Phenotype Files** (`--phenotypes_dir`): Directory containing TSV files, each with format:
   ```
   ID      Phenotype_Value
   Sample001       5.2
   Sample002       3.8
   Sample003       -9
   ```
   - Must have exactly 2 columns
   - First column must be named "ID"
   - Missing values should be "-9"
   - One file per phenotype

3. **VCF Files** (at least one of):
   - `--snp_vcf`: SNP variants
   - `--indel_vcf`: INDEL variants
   - `--sv_vcf`: Structural variants

4. **Output Configuration**:
   - `--output_prefix`: Prefix for output directories (e.g., "532rice")
   - `--analysis_type`: Analysis type name (e.g., "SNP_INDEL_SV_split")

### Input Directory Structure

```
GraphPanGWAS/
├── raw_vcf_dir/
│   ├── 532rice.snp.impute.vcf
│   ├── 532rice.indel.impute.vcf
│   └── 532rice.pangenie.multi.impute.vcf
├── phenotypes/
│   └── 532rice/
│       ├── grain_length.tsv
│       ├── grain_width.tsv
│       └── plant_height.tsv
└── 532rice.list
```

## Usage

### Quick Start

```bash
# Local execution
nextflow run main.nf \
  --population_list 532rice.list \
  --snp_vcf raw_vcf_dir/532rice.snp.impute.vcf \
  --indel_vcf raw_vcf_dir/532rice.indel.impute.vcf \
  --sv_vcf raw_vcf_dir/532rice.pangenie.multi.impute.vcf \
  --phenotypes_dir ./phenotypes/532rice \
  --output_prefix 532rice \
  --analysis_type SNP_INDEL_SV_split \
  -profile local

# SLURM cluster execution
nextflow run main.nf \
  --population_list 532rice.list \
  --snp_vcf raw_vcf_dir/532rice.snp.impute.vcf \
  --indel_vcf raw_vcf_dir/532rice.indel.impute.vcf \
  --sv_vcf raw_vcf_dir/532rice.pangenie.multi.impute.vcf \
  --phenotypes_dir ./phenotypes/532rice \
  --output_prefix 532rice \
  --analysis_type SNP_INDEL_SV_split \
  -profile slurm
```

### Resume Failed Runs

Nextflow supports automatic resume of failed pipelines:

```bash
nextflow run main.nf \
  --population_list 532rice.list \
  ... [other parameters] \
  -profile slurm \
  -resume
```

### Advanced Options

```bash
nextflow run main.nf \
  --population_list 532rice.list \
  --snp_vcf raw_vcf_dir/532rice.snp.impute.vcf \
  --phenotypes_dir ./phenotypes/532rice \
  --output_prefix 532rice \
  --analysis_type SNP_only \
  --max_parallel_gwas 10 \
  --maf_threshold 0.05 \
  --clump_p1 1e-5 \
  --clump_r2 0.2 \
  -profile slurm
```

## Configuration

### Execution Profiles

The pipeline includes two execution profiles:

#### Local Profile (`-profile local`)
- Runs all processes on the local machine
- Suitable for small datasets or testing
- Limited parallelization (max 5 concurrent GWAS jobs)

#### SLURM Profile (`-profile slurm`)
- Submits jobs to SLURM cluster
- Automatic resource allocation and queuing
- High parallelization (up to 50 concurrent jobs)
- Edit `conf/slurm.config` to customize queue names and resource limits

### Customizing Tool Paths

If tools are not in your PATH, specify their locations:

```bash
nextflow run main.nf \
  --gemma /path/to/gemma \
  --plink /path/to/plink \
  --gcta /path/to/gcta64 \
  ... [other parameters]
```

Or edit `conf/params.config`:

```groovy
params {
    gemma    = '/path/to/gemma'
    plink    = '/path/to/plink'
    gcta     = '/path/to/gcta64'
    bcftools = '/path/to/bcftools'
}
```

### Analysis Parameters

Key parameters can be adjusted in `conf/params.config` or via command line:

```bash
--maf_threshold 0.01          # Minor allele frequency threshold
--miss_threshold 0.2          # Missing rate threshold
--clump_p1 1e-6              # Clumping p-value threshold
--clump_r2 0.1               # Clumping LD r² threshold
--clump_kb 1000              # Clumping window size (kb)
--max_parallel_gwas 5        # Max parallel GWAS jobs
```

## Output Structure

```
532rice/
├── SNP_INDEL_SV_split/
│   ├── results/
│   │   ├── grain_length/
│   │   │   ├── gemma_lmm.assoc.txt    # Association results
│   │   │   └── grain_length.hsq        # Heritability estimate
│   │   ├── grain_width/
│   │   └── plant_height/
│   ├── lead_markers/
│   │   ├── grain_length.clumped        # Independent lead SNPs
│   │   ├── grain_width.clumped
│   │   └── plant_height.clumped
│   ├── heritability/
│   │   └── 532rice.SNP_INDEL_SV_split.gcta.tsv
│   ├── Lead_marker/
│   │   └── 532rice.SNP_INDEL_SV_split.ld_marker.tsv
│   └── GWAS_summary_report.txt
├── pipeline_info/
│   ├── execution_timeline.html         # Execution timeline
│   ├── execution_report.html           # Resource usage report
│   ├── execution_trace.txt             # Detailed trace
│   └── pipeline_dag.svg                # Pipeline DAG visualization
└── vcf_dir/
    └── merged.SNP_INDEL_SV_split.vcf   # Processed merged VCF
```

## Key Output Files

1. **Association Results** (`results/[phenotype]/gemma_lmm.assoc.txt`):
   - Genome-wide association test results
   - Columns: chr, rs, ps, n_miss, allele1, allele0, af, beta, se, logl_H1, l_remle, l_mle, p_wald, p_lrt, p_score

2. **Heritability Estimates** (`results/[phenotype]/[phenotype].hsq`):
   - GCTA-GREML heritability analysis
   - V(G)/Vp: SNP-based heritability

3. **Lead Markers** (`lead_markers/[phenotype].clumped`):
   - Independent associated variants after LD clumping
   - Index SNPs representing distinct association signals

4. **Summary Files**:
   - `heritability/[prefix].[type].gcta.tsv`: All heritability estimates
   - `Lead_marker/[prefix].[type].ld_marker.tsv`: All lead markers

## Troubleshooting

### Common Issues

1. **"Command not found" errors**
   - Ensure all required tools are installed and in PATH
   - Or specify full paths using `--gemma`, `--plink`, etc.

2. **"File not found" errors**
   - Check that VCF files are in `raw_vcf_dir/`
   - Verify phenotype files are in the specified `--phenotypes_dir`
   - Ensure population list file exists

3. **Memory errors**
   - Adjust memory limits in `conf/local.config` or `conf/slurm.config`
   - Increase `--plink_memory` parameter

4. **SLURM queue issues**
   - Edit `conf/slurm.config` to match your cluster's queue names
   - Adjust `clusterOptions` if needed

### Checking Pipeline Status

Monitor your pipeline:

```bash
# View Nextflow log
tail -f .nextflow.log

# List work directories
ls -lh work/

# Check SLURM jobs (if using SLURM)
squeue -u $USER
```

### Cleaning Up

Remove intermediate files:

```bash
# Remove work directory (keeps results)
rm -rf work/

# Remove all Nextflow metadata
nextflow clean -f
```

## Performance Tips

1. **Adjust parallel GWAS jobs**: Set `--max_parallel_gwas` based on your resources
   - Local: 2-5 jobs
   - SLURM: 10-50 jobs

2. **Use SLURM for large datasets**: Much faster than local execution

3. **Resume failed runs**: Always use `-resume` to avoid recomputing completed tasks

4. **Monitor resources**: Check `pipeline_info/execution_report.html` to optimize resource allocation

## Comparison with Original Bash Pipeline

### Advantages of Nextflow Version

✅ **Automatic parallelization**: Run multiple phenotype analyses simultaneously  
✅ **Resume capability**: Restart from failure point without re-running completed tasks  
✅ **Resource management**: Automatic SLURM job submission and monitoring  
✅ **Reproducibility**: Track all parameters, software versions, and execution details  
✅ **Scalability**: Easily scale from local to HPC cluster  
✅ **Error handling**: Automatic retries and better error reporting  
✅ **Portability**: Run on any system with Nextflow (local, SLURM, PBS, AWS, etc.)  

### Migration from Bash Scripts

The original bash scripts (`main.sh` and `gwas.sh`) are preserved for reference. Key differences:

1. **VCF preprocessing**: Now modular processes in `modules/vcf_preprocessing.nf`
2. **GWAS execution**: Automated parallel execution in `modules/gwas.nf`
3. **Configuration**: Centralized in `conf/*.config` files
4. **Logging**: Automatic execution reports and resource usage tracking

## Citation

If you use this pipeline, please cite:

- **GEMMA**: Zhou X, Stephens M (2012) Nature Genetics 44:821-824
- **PLINK**: Purcell S, et al. (2007) Am J Hum Genet 81:559-575
- **GCTA**: Yang J, et al. (2011) Am J Hum Genet 88:76-82
- **Nextflow**: Di Tommaso P, et al. (2017) Nature Biotechnology 35:316-319

## Support

For issues and questions:
- Check the troubleshooting section above
- Review execution reports in `pipeline_info/`
- Examine work directory logs for specific failures

## License

[Add your license information here]

