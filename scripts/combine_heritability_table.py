#!/usr/bin/env python3
"""
Combine multiple LDAK heritability TSV files into a single Excel table.

This script reads all ${output_prefix}.*.ldak.tsv files and combines them into
a formatted Excel table with hierarchical column headers:
- Sheet 1 (Heritability): Summary table with all phenotypes
  - Phenotype (first column)
  - Heritability | split | SNP/INDEL/SV/SNP_INDEL/SNP_INDEL_SV
  - Heritability | Unsplit | SNP/INDEL/SV/SNP_INDEL/SNP_INDEL_SV
- Additional sheets for merged types (SNP_INDEL and SNP_INDEL_SV):
  - Detailed REML information including Converged status
  - Her_K1, Her_K2, Her_K3, Her_Top, Her_All with SE values
"""

import sys
import os
import glob
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill

def parse_way(filename, output_prefix):
    """Extract way type from filename like 'prefix.SNP_split.ldak.tsv'"""
    basename = os.path.basename(filename)
    # Remove output_prefix and .ldak.tsv
    way_part = basename.replace(f"{output_prefix}.", "").replace(".ldak.tsv", "")
    return way_part

def read_tsv_file(filepath):
    """Read a single TSV file and return DataFrame"""
    try:
        df = pd.read_csv(filepath, sep='\t')
        return df
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}", file=sys.stderr)
        return None

def parse_reml_file(filepath):
    """Parse LDAK REML file and extract key information"""
    try:
        data = {}
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Parse header information
        for line in lines:
            line = line.strip()
            if line.startswith('Converged'):
                parts = line.split()
                data['Converged'] = parts[1] if len(parts) > 1 else 'NA'
            elif line.startswith('Component'):
                # This is the header line for component data
                continue
            elif line.startswith('Her_K1'):
                parts = line.split()
                if len(parts) >= 3:
                    data['Her_K1_Heritability'] = float(parts[1]) if parts[1] != 'NA' else None
                    data['Her_K1_SE'] = float(parts[2]) if parts[2] != 'NA' else None
            elif line.startswith('Her_K2'):
                parts = line.split()
                if len(parts) >= 3:
                    data['Her_K2_Heritability'] = float(parts[1]) if parts[1] != 'NA' else None
                    data['Her_K2_SE'] = float(parts[2]) if parts[2] != 'NA' else None
            elif line.startswith('Her_K3'):
                parts = line.split()
                if len(parts) >= 3:
                    data['Her_K3_Heritability'] = float(parts[1]) if parts[1] != 'NA' else None
                    data['Her_K3_SE'] = float(parts[2]) if parts[2] != 'NA' else None
            elif line.startswith('Her_Top'):
                parts = line.split()
                if len(parts) >= 3:
                    data['Her_Top_Heritability'] = float(parts[1]) if parts[1] != 'NA' else None
                    data['Her_Top_SE'] = float(parts[2]) if parts[2] != 'NA' else None
            elif line.startswith('Her_All'):
                parts = line.split()
                if len(parts) >= 3:
                    data['Her_All_Heritability'] = float(parts[1]) if parts[1] != 'NA' else None
                    data['Her_All_SE'] = float(parts[2]) if parts[2] != 'NA' else None
        
        return data
    except Exception as e:
        print(f"Warning: Could not parse REML file {filepath}: {e}", file=sys.stderr)
        return None

def find_reml_files(directory, way_type):
    """Find all REML files for a specific way type (e.g., SNP_INDEL_split or SNP_INDEL_SV_unsplit)"""
    reml_files = {}
    
    # Pattern: directory/way_type/results/phenotype/*.reml
    pattern = os.path.join(directory, way_type, "results", "*", "*.reml")
    files = glob.glob(pattern)
    
    for filepath in files:
        # Extract phenotype from path: .../results/phenotype_name/file.reml
        phenotype = os.path.basename(os.path.dirname(filepath))
        reml_files[phenotype] = filepath
    
    return reml_files

def add_reml_detail_sheet(wb, directory, way_type, sheet_name, phenotypes):
    """Add a sheet with detailed REML information for a specific way type"""
    
    # Find all REML files for this way type
    reml_files = find_reml_files(directory, way_type)
    
    if not reml_files:
        print(f"Warning: No REML files found for {way_type}", file=sys.stderr)
        return
    
    # Create new sheet
    ws = wb.create_sheet(title=sheet_name)
    
    # Define headers
    headers = ['Phenotype', 'Converged', 
               'Her_K1', 'Her_K1_SE',
               'Her_K2', 'Her_K2_SE',
               'Her_K3', 'Her_K3_SE',
               'Her_Top', 'Her_Top_SE',
               'Her_All', 'Her_All_SE']
    
    # Write headers
    for col_idx, header in enumerate(headers, start=1):
        cell = ws.cell(row=1, column=col_idx, value=header)
        cell.font = Font(bold=True, size=11)
        cell.fill = PatternFill(start_color='D3D3D3', end_color='D3D3D3', fill_type='solid')
        cell.alignment = Alignment(horizontal='center', vertical='center')
    
    # Write data for each phenotype
    row_idx = 2
    for phenotype in phenotypes:
        if phenotype in reml_files:
            reml_data = parse_reml_file(reml_files[phenotype])
            if reml_data:
                ws.cell(row=row_idx, column=1, value=phenotype)
                ws.cell(row=row_idx, column=2, value=reml_data.get('Converged', 'NA'))
                
                # Her_K1
                ws.cell(row=row_idx, column=3, value=reml_data.get('Her_K1_Heritability'))
                ws.cell(row=row_idx, column=4, value=reml_data.get('Her_K1_SE'))
                
                # Her_K2
                ws.cell(row=row_idx, column=5, value=reml_data.get('Her_K2_Heritability'))
                ws.cell(row=row_idx, column=6, value=reml_data.get('Her_K2_SE'))
                
                # Her_K3
                ws.cell(row=row_idx, column=7, value=reml_data.get('Her_K3_Heritability'))
                ws.cell(row=row_idx, column=8, value=reml_data.get('Her_K3_SE'))
                
                # Her_Top
                ws.cell(row=row_idx, column=9, value=reml_data.get('Her_Top_Heritability'))
                ws.cell(row=row_idx, column=10, value=reml_data.get('Her_Top_SE'))
                
                # Her_All
                ws.cell(row=row_idx, column=11, value=reml_data.get('Her_All_Heritability'))
                ws.cell(row=row_idx, column=12, value=reml_data.get('Her_All_SE'))
                
                # Format numeric cells
                for col in range(3, 13):
                    cell = ws.cell(row=row_idx, column=col)
                    if cell.value is not None and isinstance(cell.value, (int, float)):
                        cell.number_format = '0.000000'
                
                row_idx += 1
    
    # Adjust column widths
    ws.column_dimensions['A'].width = 30
    ws.column_dimensions['B'].width = 12
    for col_letter in ['C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']:
        ws.column_dimensions[col_letter].width = 14
    
    print(f"Added sheet '{sheet_name}' with {row_idx - 2} phenotypes", file=sys.stderr)

def combine_heritability_files(directory, output_prefix, output_file):
    """Combine all heritability TSV files into Excel table"""
    
    # Find all matching TSV files
    pattern = os.path.join(directory, f"{output_prefix}.*.ldak.tsv")
    tsv_files = glob.glob(pattern)
    
    if not tsv_files:
        print(f"Error: No files found matching pattern {pattern}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(tsv_files)} TSV files to combine", file=sys.stderr)
    
    # Read all TSV files and store by way type
    data_dict = {}  # way -> DataFrame
    all_phenotypes = set()
    
    for tsv_file in tsv_files:
        way = parse_way(tsv_file, output_prefix)
        df = read_tsv_file(tsv_file)
        if df is not None and not df.empty:
            data_dict[way] = df
            all_phenotypes.update(df['Phenotype'].tolist())
    
    if not data_dict:
        print("Error: No valid data found in any TSV file", file=sys.stderr)
        sys.exit(1)
    
    # Define the order of way types (matching the image)
    way_types = ['SNP', 'INDEL', 'SV', 'SNP_INDEL', 'SNP_INDEL_SV']
    split_types = [f"{wt}_split" for wt in way_types]
    unsplit_types = [f"{wt}_unsplit" for wt in way_types]
    
    # Create ordered list of phenotypes
    sorted_phenotypes = sorted(all_phenotypes)
    
    # Build data structure for Excel output
    # We'll create a structure that matches the hierarchical header format
    data_rows = []
    
    # Process each phenotype
    for phenotype in sorted_phenotypes:
        row = [phenotype]  # Start with Phenotype column
        
        # Process split types: only Heritability (no P_value)
        for way in split_types:
            if way in data_dict:
                df = data_dict[way]
                pheno_data = df[df['Phenotype'] == phenotype]
                if not pheno_data.empty:
                    row.append(pheno_data.iloc[0]['Heritability'])
                else:
                    row.append(None)
            else:
                row.append(None)
        
        # Process unsplit types: only Heritability (no P_value)
        for way in unsplit_types:
            if way in data_dict:
                df = data_dict[way]
                pheno_data = df[df['Phenotype'] == phenotype]
                if not pheno_data.empty:
                    row.append(pheno_data.iloc[0]['Heritability'])
                else:
                    row.append(None)
            else:
                row.append(None)
        
        data_rows.append(row)
    
    # Calculate mean row
    mean_row = ['mean']
    for col_idx in range(1, len(data_rows[0])):  # Skip Phenotype column
        col_values = [row[col_idx] for row in data_rows if row[col_idx] is not None and pd.notna(row[col_idx])]
        if col_values:
            mean_val = sum(col_values) / len(col_values)
            mean_row.append(mean_val)
        else:
            mean_row.append(None)
    
    data_rows.append(mean_row)
    
    # Create Excel file with hierarchical headers
    wb = Workbook()
    ws = wb.active
    ws.title = "Heritability"
    
    # Write headers (2 rows) - no P_value columns
    # Row 1: Main header
    ws.cell(row=1, column=1, value='Phenotype')
    
    # Calculate positions for split section (only Heritability, no P_value)
    split_start_col = 2
    split_end_col = split_start_col + len(split_types) - 1
    
    ws.cell(row=1, column=split_start_col, value='Heritability')
    ws.merge_cells(start_row=1, start_column=split_start_col, end_row=1, end_column=split_end_col)
    
    unsplit_start_col = split_end_col + 1
    unsplit_end_col = unsplit_start_col + len(unsplit_types) - 1
    
    ws.cell(row=1, column=unsplit_start_col, value='Heritability')
    ws.merge_cells(start_row=1, start_column=unsplit_start_col, end_row=1, end_column=unsplit_end_col)
    
    # Row 2: split/Unsplit sub-headers
    ws.cell(row=2, column=split_start_col, value='split')
    ws.merge_cells(start_row=2, start_column=split_start_col, end_row=2, end_column=split_end_col)
    
    ws.cell(row=2, column=unsplit_start_col, value='Unsplit')
    ws.merge_cells(start_row=2, start_column=unsplit_start_col, end_row=2, end_column=unsplit_end_col)
    
    # Row 3: Individual way types (one column per type, no P_value)
    col = split_start_col
    for way_type in way_types:
        ws.cell(row=3, column=col, value=way_type)
        col += 1
    
    col = unsplit_start_col
    for way_type in way_types:
        ws.cell(row=3, column=col, value=way_type)
        col += 1
    
    # Write data rows (starting from row 4)
    for row_idx, data_row in enumerate(data_rows, start=4):
        for col_idx, value in enumerate(data_row, start=1):
            cell = ws.cell(row=row_idx, column=col_idx, value=value)
            # Format numeric cells - all are Heritability values now
            if isinstance(value, (int, float)) and pd.notna(value) and col_idx > 1:
                cell.number_format = '0.000000'
    
    # Format headers
    header_fill = PatternFill(start_color='D3D3D3', end_color='D3D3D3', fill_type='solid')
    header_font = Font(bold=True, size=11)
    center_alignment = Alignment(horizontal='center', vertical='center')
    
    for row in range(1, 4):
        for col in range(1, ws.max_column + 1):
            cell = ws.cell(row=row, column=col)
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = center_alignment
    
    # Merge Phenotype header across first 3 rows
    ws.merge_cells(start_row=1, start_column=1, end_row=3, end_column=1)
    
    # Adjust column widths
    ws.column_dimensions['A'].width = 30
    for col in range(2, ws.max_column + 1):
        ws.column_dimensions[chr(64 + col)].width = 12
    
    # Add detailed REML sheets for merged types (SNP_INDEL and SNP_INDEL_SV)
    merged_types = [
        ('SNP_INDEL_split', 'SNP_INDEL_split'),
        ('SNP_INDEL_unsplit', 'SNP_INDEL_unsplit'),
        ('SNP_INDEL_SV_split', 'SNP_INDEL_SV_split'),
        ('SNP_INDEL_SV_unsplit', 'SNP_INDEL_SV_unsplit')
    ]
    
    for way_type, sheet_name in merged_types:
        add_reml_detail_sheet(wb, directory, way_type, sheet_name, sorted_phenotypes)
    
    # Save file
    wb.save(output_file)
    
    print(f"Successfully created Excel file: {output_file}", file=sys.stderr)

def main():
    if len(sys.argv) != 4:
        print("Usage: combine_heritability_table.py <directory> <output_prefix> <output_file>", file=sys.stderr)
        sys.exit(1)
    
    directory = sys.argv[1]
    output_prefix = sys.argv[2]
    output_file = sys.argv[3]
    
    if not os.path.isdir(directory):
        print(f"Error: Directory {directory} does not exist", file=sys.stderr)
        sys.exit(1)
    
    combine_heritability_files(directory, output_prefix, output_file)

if __name__ == '__main__':
    main()

