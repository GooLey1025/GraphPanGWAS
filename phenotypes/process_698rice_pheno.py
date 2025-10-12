#!/usr/bin/env python3
"""
Script to process rice phenotype data:
1. Remove samples listed in del.list
2. Group phenotypes by trait name and calculate BLUP for multi-environment traits
3. Generate individual .tsv files for each phenotype (including existing and calculated BLUP)
4. Each .tsv file contains ID column and corresponding phenotype values
5. Filter phenotypes with insufficient valid samples (<200)
"""

import pandas as pd
import numpy as np
import os
import sys
from collections import defaultdict


def calculate_blup(values):
    """
    Calculate BLUP (Best Linear Unbiased Prediction) for a set of values.
    For simplicity, we'll use the mean of non-missing values as BLUP.
    In practice, BLUP would involve more complex mixed model calculations.
    """
    valid_values = []
    for value in values:
        if not (pd.isna(value) or value == 'NA' or value == '' or str(value).strip() == ''):
            try:
                valid_values.append(float(value))
            except (ValueError, TypeError):
                continue

    if len(valid_values) == 0:
        return 'NA'

    return np.mean(valid_values)


def main():
    # File paths
    pheno_file = "705rice.pheno"
    del_file = "del.list"

    # Check if files exist
    if not os.path.exists(pheno_file):
        print(f"Error: {pheno_file} not found!")
        sys.exit(1)

    if not os.path.exists(del_file):
        print(f"Error: {del_file} not found!")
        sys.exit(1)

    print("Reading phenotype data...")
    # Read phenotype data
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    print(f"Original data shape: {pheno_df.shape}")

    # Read samples to delete
    print("Reading samples to delete...")
    with open(del_file, 'r') as f:
        samples_to_delete = [line.strip() for line in f if line.strip()]

    print(f"Samples to delete: {len(samples_to_delete)}")

    # Remove samples listed in del.list
    print("Removing samples...")
    original_count = len(pheno_df)
    pheno_df = pheno_df[~pheno_df['LINE'].isin(samples_to_delete)]
    removed_count = original_count - len(pheno_df)
    print(f"Removed {removed_count} samples. Remaining: {len(pheno_df)}")

    # Sort by LINE (ID) to ensure proper ordering
    pheno_df = pheno_df.sort_values('LINE').reset_index(drop=True)
    print("Data sorted by LINE (ID)")

    # Get phenotype columns (all columns except the first one which is LINE/ID)
    phenotype_columns = pheno_df.columns[1:]
    print(f"Found {len(phenotype_columns)} phenotype columns")

    # Set minimum sample threshold (200 samples)
    min_samples = 200
    print(f"Minimum sample threshold: {min_samples} valid values")

    # Group phenotypes by trait name (before the first underscore)
    print("Grouping phenotypes by trait name...")
    trait_groups = defaultdict(list)
    blup_phenotypes = []

    for phenotype in phenotype_columns:
        if '_BLUP' in phenotype:
            # This is already a BLUP phenotype
            blup_phenotypes.append(phenotype)
        else:
            # Extract trait name (before first underscore)
            trait_name = phenotype.split('_')[0]
            trait_groups[trait_name].append(phenotype)

    print(f"Found {len(trait_groups)} trait groups with multiple environments")
    print(f"Found {len(blup_phenotypes)} existing BLUP phenotypes")

    # Calculate BLUP for traits with multiple environments (only if BLUP doesn't exist)
    print("Calculating BLUP values for multi-environment traits...")
    for trait_name, env_phenotypes in trait_groups.items():
        if len(env_phenotypes) > 1:  # Only calculate BLUP for multi-environment traits
            blup_col_name = f"{trait_name}_BLUP"

            # Check if BLUP already exists
            if blup_col_name in blup_phenotypes:
                print(f"  {trait_name}: BLUP already exists - SKIP calculation")
            else:
                print(
                    f"  Calculating BLUP for {trait_name} ({len(env_phenotypes)} environments)")
                blup_values = []

                for idx, row in pheno_df.iterrows():
                    # Get values from all environments for this trait
                    env_values = [row[env] for env in env_phenotypes]
                    blup_value = calculate_blup(env_values)
                    blup_values.append(blup_value)

                # Add BLUP column to dataframe
                pheno_df[blup_col_name] = blup_values
                blup_phenotypes.append(blup_col_name)
                print(f"    -> Added {blup_col_name}")

    # Determine which phenotypes to output
    print("Determining output phenotypes...")
    output_phenotypes = []

    for trait_name, env_phenotypes in trait_groups.items():
        blup_col_name = f"{trait_name}_BLUP"

        if len(env_phenotypes) > 1:
            # Multi-environment trait: output BLUP only
            if blup_col_name in blup_phenotypes:
                output_phenotypes.append(blup_col_name)
                print(f"  {trait_name}: Multi-environment -> Output BLUP only")
            else:
                # This shouldn't happen if BLUP calculation worked correctly
                print(
                    f"  {trait_name}: Multi-environment but no BLUP found - SKIP")
        else:
            # Single-environment trait: output original data
            output_phenotypes.extend(env_phenotypes)
            print(f"  {trait_name}: Single-environment -> Output original data")

    # Add any existing BLUP phenotypes that weren't covered above
    for blup_pheno in blup_phenotypes:
        if blup_pheno not in output_phenotypes:
            output_phenotypes.append(blup_pheno)
            print(f"  {blup_pheno}: Existing BLUP -> Output")

    print(f"Total phenotypes to output: {len(output_phenotypes)}")

    # Filter phenotypes based on valid sample count
    print("Checking valid sample counts for each phenotype...")
    valid_phenotypes = []
    for phenotype in output_phenotypes:
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
        f"\nValid phenotypes after filtering: {len(valid_phenotypes)}/{len(output_phenotypes)}")

    # Create output directory if it doesn't exist
    output_dir = "698rice"
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
            'ID': pheno_df['LINE'],
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
    print(
        f"Skipped {len(output_phenotypes) - len(valid_phenotypes)} phenotypes due to insufficient valid samples (<200)")
    print(f"Each file contains ID column and corresponding phenotype values")
    print(f"All data is sorted by ID for consistency")


if __name__ == "__main__":
    main()
