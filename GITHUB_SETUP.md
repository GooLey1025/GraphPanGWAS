# Setting Up GitHub Repository for GraphPanGWAS

## Step-by-Step Instructions

### 1. Initialize Git Repository Locally

```bash
cd /data/xhhuanglab/gulei/projects/GraphPanGWAS

# Initialize git
git init

# Add all files
git add .

# Make initial commit
git commit -m "Initial commit: Nextflow GWAS pipeline with split/unsplit analysis"
```

### 2. Create GitHub Repository

#### Option A: Using GitHub Web Interface (Easier)

1. Go to https://github.com
2. Click the "+" icon in top right → "New repository"
3. Repository name: `GraphPanGWAS`
4. Description: `Nextflow pipeline for GWAS analysis of pangenome variants (SNP, INDEL, SV)`
5. Choose: **Public** or **Private**
6. **DO NOT** initialize with README (we already have one)
7. Click "Create repository"

#### Option B: Using GitHub CLI (if installed)

```bash
gh repo create GraphPanGWAS --public --source=. --remote=origin --push
```

### 3. Connect Local Repository to GitHub

After creating the repository on GitHub, connect it:

```bash
# Add GitHub remote (replace YOUR_USERNAME with your GitHub username)
git remote add origin https://github.com/YOUR_USERNAME/GraphPanGWAS.git

# Verify remote
git remote -v

# Push to GitHub
git branch -M main
git push -u origin main
```

### 4. Verify Upload

Go to: `https://github.com/YOUR_USERNAME/GraphPanGWAS`

You should see:
- All Nextflow files (main.nf, modules/, conf/)
- Documentation (README.md, QUICKSTART.md, etc.)
- Configuration files (params.yaml, nextflow.config)
- Scripts directory

### 5. Update Repository (Future Changes)

```bash
# Stage changes
git add .

# Commit with message
git commit -m "Description of changes"

# Push to GitHub
git push
```

## What Gets Uploaded

### ✅ Included (Pipeline Code):
- `main.nf`
- `nextflow.config`
- `params.yaml`
- `modules/*.nf`
- `conf/*.config`
- `scripts/*.sh`, `scripts/*.py`, `scripts/*.R`
- Documentation: `*.md` files
- Utility scripts: `validate_setup.sh`, `run_example.sh`

### ❌ Excluded (Data & Results):
- `work/` - Nextflow work directory
- `*.vcf`, `*.vcf.gz` - VCF files (too large)
- `phenotypes/` - Phenotype data (optional, can be included)
- Output directories (results, heritability, etc.)
- Log files

## GitHub Best Practices

### Add a License

```bash
# Create LICENSE file (MIT License example)
cat > LICENSE << 'EOFLIC'
MIT License

Copyright (c) 2025 [Your Name]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOFLIC

git add LICENSE
git commit -m "Add MIT license"
git push
```

### Add Topics/Tags on GitHub

After creating the repository, add topics on GitHub:
1. Go to your repository page
2. Click "⚙️" next to "About"
3. Add topics: `nextflow`, `gwas`, `bioinformatics`, `genomics`, `pipeline`, `ngs`

### Enable GitHub Actions (Optional)

Create `.github/workflows/test.yml` for automated testing.

## Troubleshooting

### Issue: Permission denied (publickey)

Use HTTPS instead of SSH:
```bash
git remote set-url origin https://github.com/YOUR_USERNAME/GraphPanGWAS.git
```

### Issue: Files too large

GitHub has 100MB file limit. Check with:
```bash
find . -size +50M -not -path "./.git/*"
```

Add large files to `.gitignore`

### Issue: Authentication failed

Use personal access token:
1. GitHub → Settings → Developer settings → Personal access tokens
2. Generate new token with `repo` scope
3. Use token as password when pushing

## Example Complete Setup

```bash
# Navigate to project
cd /data/xhhuanglab/gulei/projects/GraphPanGWAS

# Initialize and commit
git init
git add .
git commit -m "Initial commit: GraphPanGWAS Nextflow pipeline"

# Connect to GitHub (replace YOUR_USERNAME)
git remote add origin https://github.com/YOUR_USERNAME/GraphPanGWAS.git
git branch -M main
git push -u origin main

# Done! 🎉
```

## After GitHub Setup

Share your repository:
```
https://github.com/YOUR_USERNAME/GraphPanGWAS
```

Users can clone and use:
```bash
git clone https://github.com/YOUR_USERNAME/GraphPanGWAS.git
cd GraphPanGWAS
nextflow run main.nf -params-file params.yaml -profile slurm
```
