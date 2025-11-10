#!/usr/bin/env python3
"""
Manhattan Plot Generator for GWAS Results
Generates a Manhattan plot from GEMMA association results
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib import font_manager
from matplotlib.lines import Line2D

# Set matplotlib parameters for better looking plots
font_dir = os.path.expanduser("~/.local/share/fonts/arial")
if os.path.isdir(font_dir):
    for f in font_manager.findSystemFonts(fontpaths=[font_dir], fontext='ttf'):
        font_manager.fontManager.addfont(f)

if any('arial' in font.name.lower() for font in font_manager.fontManager.ttflist):
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
else:
    rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10
rcParams['axes.linewidth'] = 1.5


def preprocess_label(raw_label):
    """
    Clean phenotype label for axis display.
    """
    if not raw_label:
        return ""

    label = raw_label.replace('processed_', '')
    suffixes = [
        '_SNP_split',
        '_INDEL_split',
        '_SNP_INDEL_SV_split',
        '_SNP_INDEL_SV_unsplit',
    ]
    for suffix in suffixes:
        if label.endswith(suffix):
            label = label[: -len(suffix)]
            break

    return label


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

    # Derive variant type from the second column (e.g., marker ID)
    second_col = df.columns[1]
    variant_values = df[second_col].astype(str)
    df['variant_type'] = 'SNP'
    indel_mask = variant_values.str.contains('INDEL', case=False, na=False)
    sv_mask = variant_values.str.contains('SV', case=False, na=False)
    df.loc[indel_mask, 'variant_type'] = 'INDEL'
    df.loc[sv_mask, 'variant_type'] = 'SV'
    
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


def create_manhattan_plot(df, output_file, label, significance_threshold):
    """
    Create Manhattan plot from GWAS data
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with columns: chr, ps, log_p
    output_file : str
        Path to output PNG file
    significance_threshold : float
        P-value threshold for genome-wide significance
    label : str
        Label for the plot
    """
    
    # Calculate cumulative position for each marker
    df['cumulative_pos'] = 0
    cumulative_offset = 0
    chr_centers = []
    chr_boundaries = [0]
    
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
        chr_boundaries.append(cumulative_offset)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 2.5))
    
    # Plot points by variant type
    variant_styles = {
        'SNP': {'color': '#1f77b4', 'marker': 'o', 'size': 9},
        'INDEL': {'color': '#2ca02c', 'marker': '^', 'size': 20},
        'SV': {'color': '#ff7f0e', 'marker': '*', 'size': 36},
    }
    
    legend_handles = []
    
    for variant_type, style in variant_styles.items():
        subset = df[df['variant_type'] == variant_type]
        if subset.empty:
            continue
        
        ax.scatter(
            subset['cumulative_pos'],
            subset['log_p'],
            c=style['color'],
            s=style['size'],
            marker=style['marker'],
            alpha=0.85,
            linewidths=0,
        )
        
        legend_handles.append(
            Line2D(
                [],
                [],
                linestyle='',
                marker=style['marker'],
                markerfacecolor=style['color'],
                markeredgecolor=style['color'],
                markersize=6,
                label=variant_type,
            )
        )
    
    # Add genome-wide significance line
    sig_line = -np.log10(significance_threshold) if significance_threshold > 0 else 0
    ax.axhline(y=sig_line, color='black', linestyle='--', linewidth=1.2, dashes=(4, 4), alpha=0.7)
    
    # Draw chromosome boundary lines
    for boundary in chr_boundaries[1:-1]:
        ax.axvline(x=boundary, color='#d3d3d3', linestyle='-', linewidth=0.8, alpha=0.7)
    
    # Labels and title
    if label != "":
        axis_label = f'Chromosome - {label}'
    else:
        axis_label = 'Chromosome'
    ax.set_xlabel(axis_label, fontsize=12, fontweight='normal')
    ax.set_ylabel(r'$-\log_{10}[\mathit{P}]$', fontsize=12, fontweight='normal')
    
    # Set y-axis limits
    y_max = max(df['log_p'].max() * 1.1, sig_line * 1.2 if significance_threshold > 0 else 0)
    y_max *= 1.05
    ax.set_ylim(0, y_max)
    
    # Set x-axis limits to remove gaps
    ax.set_xlim(df['cumulative_pos'].min(), df['cumulative_pos'].max())
    
    # Annotate chromosome numbers above intervals
    chrom_label_y = y_max * 0.92
    for chrom, center in chr_centers:
        ax.text(center, chrom_label_y, str(int(chrom)), ha='center', va='bottom', fontfamily='Arial', fontsize=12)
    
    ax.set_xticks([])
    ax.tick_params(axis='x', which='both', length=0)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add legend for variant types
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc='upper left',
            bbox_to_anchor=(0, 0.92),
            frameon=False,
            fontsize=11,
            handlelength=1.2,
        )
    
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

    cleaned_label = preprocess_label(label)
    df = load_gwas_data(input_file)
    num_markers = len(df)
    significance_threshold = 0.05 / num_markers if num_markers > 0 else 1
    create_manhattan_plot(df, output_file, cleaned_label, significance_threshold)


if __name__ == "__main__":
    main()

