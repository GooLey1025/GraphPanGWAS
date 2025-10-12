# GraphPanGWAS - Nextflow GWAS Pipeline

A comprehensive Nextflow pipeline for genome-wide association studies (GWAS) of pangenome variants including SNPs, INDELs, and structural variants (SVs).

## Features

- ✅ Automatic processing of both SPLIT and UNSPLIT variant analyses
- ✅ Parallel VCF preprocessing (SNP, INDEL, SV)
- ✅ Automated GWAS with GEMMA (kinship-corrected LMM)
- ✅ Lead marker identification with PLINK clumping
- ✅ SNP-based heritability estimation with GCTA
- ✅ Support for local and SLURM cluster execution
- ✅ Automatic resume capability
- ✅ Comprehensive execution reports

## Quick Start

```bash
# Run the complete pipeline
nextflow run main.nf \
  --population_list samples.list \
  --snp_vcf raw_vcf_dir/snps.vcf \
  --indel_vcf raw_vcf_dir/indels.vcf \
  --sv_vcf raw_vcf_dir/svs.vcf \
  --phenotypes_dir ./phenotypes/population \
  --output_prefix my_analysis \
  --analysis_type SNP_INDEL_SV \
  -profile slurm
```

## Installation

### Requirements

- Nextflow >= 22.10.0
- PLINK (1.9 or 2.0)
- GEMMA >= 0.98
- GCTA >= 1.93
- bcftools >= 1.10
- Python 3 with pandas

### Setup

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/GraphPanGWAS.git
cd GraphPanGWAS

# Validate setup
bash validate_setup.sh

# Run test
nextflow run main.nf -params-file params.yaml -profile local
```

## Documentation

- [QUICKSTART.md](QUICKSTART.md) - Quick start guide
- [README_NEXTFLOW.md](README_NEXTFLOW.md) - Complete documentation
- [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) - Migration from bash version
- [CHANGELOG.md](CHANGELOG.md) - Version history

## Usage

### Using Parameter File (Recommended)

```bash
# Edit params.yaml with your settings
cp params.yaml my_analysis.yaml
nano my_analysis.yaml

# Run pipeline
nextflow run main.nf -params-file my_analysis.yaml -profile slurm
```

### Command Line

```bash
nextflow run main.nf \
  --population_list samples.list \
  --phenotypes_dir ./phenotypes \
  --output_prefix output \
  --analysis_type SNP_INDEL_SV \
  --snp_vcf raw_vcf_dir/snps.vcf \
  -profile slurm
```

## Output Structure

```
output_prefix/
├── SNP_INDEL_SV_split/
│   ├── results/           # Association results
│   ├── lead_markers/      # Clumped lead SNPs
│   └── heritability/      # Heritability estimates
└── SNP_INDEL_SV_unsplit/
    ├── results/
    ├── lead_markers/
    └── heritability/
```

## Citation

If you use this pipeline, please cite:

- **GEMMA**: Zhou X, Stephens M (2012) Nature Genetics
- **PLINK**: Purcell S, et al. (2007) Am J Hum Genet
- **GCTA**: Yang J, et al. (2011) Am J Hum Genet
- **Nextflow**: Di Tommaso P, et al. (2017) Nature Biotechnology

## License

MIT License

## Contact

For issues and questions, please open a GitHub issue.
