#!/usr/bin/env python3
"""
Script to process 1171rice phenotype data:
1. Read CSV format phenotype data
2. Clean column names (remove units and special characters)
3. Generate individual .tsv files for each phenotype
4. Each .tsv file contains ID column and corresponding phenotype values
5. Filter phenotypes with insufficient valid samples (<200)
"""

import pandas as pd
import numpy as np
import os
import sys
import re


def clean_column_name(col_name):
    """
    Clean column names by removing units and special characters.
    Convert to a more standard format for file naming.
    """
    # Remove quotes if present
    col_name = col_name.strip('"')

    # Remove units in parentheses or after /
    col_name = re.sub(r'\([^)]*\)', '', col_name)  # Remove (units)
    col_name = re.sub(r'/[^/]*$', '', col_name)    # Remove /units at end
    # Remove Unicode units like <U+00B0>
    col_name = re.sub(r'<U\+[^>]*>', '', col_name)

    # Replace special characters with underscores
    col_name = re.sub(r'[^\w]', '_', col_name)

    # Remove multiple underscores and trim
    col_name = re.sub(r'_+', '_', col_name).strip('_')

    return col_name


def main():
    # File paths
    pheno_file = "1171rice.phenotype.csv"

    # Check if file exists
    if not os.path.exists(pheno_file):
        print(f"Error: {pheno_file} not found!")
        sys.exit(1)

    print("Reading 1171rice phenotype data...")
    # Read phenotype data
    pheno_df = pd.read_csv(pheno_file)
    print(f"Original data shape: {pheno_df.shape}")

    # Clean column names
    print("Cleaning column names...")
    original_columns = list(pheno_df.columns)
    cleaned_columns = []

    for col in original_columns:
        if col == "ID":
            cleaned_columns.append("ID")
        else:
            cleaned_col = clean_column_name(col)
            cleaned_columns.append(cleaned_col)
            if col != cleaned_col:
                print(f"  {col} -> {cleaned_col}")

    # Update column names
    pheno_df.columns = cleaned_columns

    # Sort by ID to ensure proper ordering
    pheno_df = pheno_df.sort_values('ID').reset_index(drop=True)
    print("Data sorted by ID")

    # Get phenotype columns (all columns except ID)
    phenotype_columns = [col for col in pheno_df.columns if col != 'ID']
    print(f"Found {len(phenotype_columns)} phenotype columns")

    # Set minimum sample threshold (200 samples)
    min_samples = 200
    print(f"Minimum sample threshold: {min_samples} valid values")

    # Filter phenotypes based on valid sample count
    print("Checking valid sample counts for each phenotype...")
    valid_phenotypes = []
    for phenotype in phenotype_columns:
        # Count valid (non-missing) values
        valid_count = 0
        total_count = len(pheno_df)

        for value in pheno_df[phenotype]:
            if not (pd.isna(value) or value == 'NA' or value == '' or str(value).strip() == ''):
                valid_count += 1

        missing_count = total_count - valid_count
        missing_rate = missing_count / total_count

        if valid_count >= min_samples:
            valid_phenotypes.append(phenotype)
            print(
                f"  {phenotype}: {valid_count} valid samples ({missing_rate:.1%} missing) - KEEP")
        else:
            print(
                f"  {phenotype}: {valid_count} valid samples ({missing_rate:.1%} missing) - SKIP (< {min_samples} samples)")

    print(
        f"\nValid phenotypes after filtering: {len(valid_phenotypes)}/{len(phenotype_columns)}")

    # Create output directory if it doesn't exist
    output_dir = "1171rice"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    # Generate individual .tsv files for each valid phenotype
    print("Generating individual .tsv files...")
    for i, phenotype in enumerate(valid_phenotypes, 1):
        print(f"Processing {i}/{len(valid_phenotypes)}: {phenotype}")

        # Create dataframe with ID and phenotype columns
        # Replace any missing values (NaN, empty strings, etc.) with 'NA'
        phenotype_values = pheno_df[phenotype].fillna('NA')
        phenotype_values = phenotype_values.replace('', 'NA')

        output_df = pd.DataFrame({
            'ID': pheno_df['ID'],
            phenotype: phenotype_values
        })

        # Output file path
        output_file = os.path.join(output_dir, f"{phenotype}.tsv")

        # Save to .tsv file
        output_df.to_csv(output_file, sep='\t', index=False)

        # Verify the file was created and show basic stats
        if os.path.exists(output_file):
            # Count non-NA values (excluding 'NA' strings)
            non_na_count = (output_df[phenotype] != 'NA').sum()
            print(f"  -> Saved: {output_file} ({non_na_count} non-NA values)")
        else:
            print(f"  -> Error: Failed to create {output_file}")

    print(f"\nProcessing complete!")
    print(
        f"Generated {len(valid_phenotypes)} .tsv files in '{output_dir}' directory")
    print(f"Skipped {len(phenotype_columns) - len(valid_phenotypes)} phenotypes due to insufficient valid samples (<200)")
    print(f"Each file contains ID column and corresponding phenotype values")
    print(f"All data is sorted by ID for consistency")


if __name__ == "__main__":
    main()
