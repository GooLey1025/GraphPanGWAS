#!/usr/bin/env python3
"""
Manhattan Plot Generator for GWAS Results
Generates a Manhattan plot from GEMMA association results
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib import rcParams

# Set matplotlib parameters for better looking plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = 10
rcParams['axes.linewidth'] = 1.5


def load_gwas_data(input_file):
    """
    Load GWAS results from GEMMA output format
    Expected columns: chr, rs, ps, n_miss, allele1, allele0, af, beta, se, logl_H1, l_remle, p_wald
    """
    # Read the association file
    df = pd.read_csv(input_file, sep='\t', header=0)
    
    # Check if required columns exist
    required_cols = ['chr', 'ps', 'p_wald']
    for col in required_cols:
        if col not in df.columns:
            sys.exit(1)
    
    # Filter out rows with invalid p-values
    df = df[df['p_wald'] > 0].copy()
    df = df[df['p_wald'] <= 1].copy()
    
    # Calculate -log10(p-value)
    df['log_p'] = -np.log10(df['p_wald'])
    
    # Convert chromosome to numeric if needed (handle chr prefixes)
    if df['chr'].dtype == 'object':
        df['chr'] = df['chr'].str.replace('chr', '', case=False)
        df['chr'] = df['chr'].str.replace('Chr', '', case=False)
    
    # Convert to numeric, removing any non-numeric chromosomes
    df['chr'] = pd.to_numeric(df['chr'], errors='coerce')
    df = df.dropna(subset=['chr'])
    df['chr'] = df['chr'].astype(int)
    
    # Sort by chromosome and position
    df = df.sort_values(['chr', 'ps'])
    
    return df


def create_manhattan_plot(df, output_file, label, significance_threshold=5e-8):
    """
    Create Manhattan plot from GWAS data
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with columns: chr, ps, log_p
    output_file : str
        Path to output PNG file
    significance_threshold : float
        P-value threshold for genome-wide significance (default: 5e-8)
    label : str
        Label for the plot
    """
    
    # Calculate cumulative position for each SNP
    df['cumulative_pos'] = 0
    cumulative_offset = 0
    chr_centers = []
    
    chromosomes = sorted(df['chr'].unique())
    
    for chrom in chromosomes:
        chr_df = df[df['chr'] == chrom]
        chr_length = chr_df['ps'].max()
        
        # Update cumulative positions
        df.loc[df['chr'] == chrom, 'cumulative_pos'] = df.loc[df['chr'] == chrom, 'ps'] + cumulative_offset
        
        # Store center position for chromosome label
        chr_centers.append((chrom, cumulative_offset + chr_length / 2))
        
        # Update offset for next chromosome
        cumulative_offset += chr_length
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Define colors for alternating chromosomes
    colors = ['#0099CC', '#FF6600']
    
    # Plot points for each chromosome with alternating colors
    for i, chrom in enumerate(chromosomes):
        chr_df = df[df['chr'] == chrom]
        color = colors[i % 2]
        
        ax.scatter(chr_df['cumulative_pos'], chr_df['log_p'], 
                  c=color, s=10, alpha=0.7, edgecolors='none')
    
    # Add genome-wide significance line
    sig_line = -np.log10(significance_threshold)
    ax.axhline(y=sig_line, color='red', linestyle='--', linewidth=1.5, dashes=(6, 3), alpha=0.5,
               label=f'Genome-wide significance (p={significance_threshold})')
    
    # Add suggestive significance line (p=1e-5)
    suggestive_line = -np.log10(1e-5)
    ax.axhline(y=suggestive_line, color='blue', linestyle='--', linewidth=1, dashes=(6, 3),
               alpha=0.5, label='Suggestive (p=1e-5)')
    
    # Set x-axis labels at chromosome centers
    chr_labels = [str(int(c[0])) for c in chr_centers]
    chr_positions = [c[1] for c in chr_centers]
    ax.set_xticks(chr_positions)
    ax.set_xticklabels(chr_labels)
    
    # Labels and title
    if label != "":
        ax.set_xlabel(f'Chromosome ({label})', fontsize=14, fontweight='bold')
    else:
        ax.set_xlabel('Chromosome', fontsize=14, fontweight='bold')
    ax.set_ylabel('-log₁₀(P)', fontsize=14, fontweight='bold')
    
    # Set y-axis limits
    y_max = max(df['log_p'].max() * 1.1, sig_line * 1.2)
    ax.set_ylim(0, y_max)
    
    # Set x-axis limits to remove gaps
    ax.set_xlim(df['cumulative_pos'].min(), df['cumulative_pos'].max())
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add legend
    ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
    
    # Add grid
    ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='png')
    
    plt.close()


def main():
    """
    Main function
    Usage: python manhattan_plot.py <input_assoc_file> <output_png> <label>
    """
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Get label if provided
    if len(sys.argv) == 4:
        if sys.argv[3] == "None":
            label = ""
        else:
            label = sys.argv[3]
    else:
        label = ""

    df = load_gwas_data(input_file)
    create_manhattan_plot(df, output_file, label)


if __name__ == "__main__":
    main()

