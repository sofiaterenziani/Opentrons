#!/usr/bin/env python
"""
DCPIP Titration Analysis Script
Analyzes absorbance data from a 384-well plate with DCPIP titrations at different pH values.

384-well plate layout:
- Columns 1-8: MES pH 6 with DCPIP concentrations [250, 150, 90, 54, 32.4, 19.44, 11.66, 0 µM]
- Columns 9-16: HEPES pH 7 with same concentrations
- Columns 17-24: HEPES pH 8 with same concentrations
- Rows: A-P (16 rows total)
"""

import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DCPIP_CONCENTRATIONS = [250, 150, 90, 54, 32.4, 19.44, 11.66, 0]  # µM

# pH conditions: pH_name -> (column_start, column_end, color)
PH_CONDITIONS = {
    'MES pH 6': (1, 8, 'blue'),
    'HEPES pH 7': (9, 16, 'green'),
    'HEPES pH 8': (17, 24, 'red')
}


def parse_absorbance_value(value):
    """Return a numeric absorbance value or NaN for overflow/missing cells."""
    if pd.isna(value):
        return np.nan

    if isinstance(value, str):
        normalized = value.strip().lower()
        if not normalized:
            return np.nan
        if normalized == 'ovrflw' or 'overflow' in normalized:
            return np.nan

    try:
        return float(value)
    except (TypeError, ValueError):
        return np.nan

def well_to_coordinates(well_label):
    """Convert well label (e.g., 'A1') to row and column indices."""
    row = ord(well_label[0]) - ord('A')  # A=0, B=1, ..., P=15
    col = int(well_label[1:]) - 1  # 1-indexed to 0-indexed
    return row, col

def parse_plate_data(df):
    """
    Parse absorbance data from the plate reader output.
    Handles multiple possible Excel formats, including the BioTek/Synergy export format used here.
    """
    data = {ph: {'conc': [], 'abs': []} for ph in PH_CONDITIONS}

    df = df.copy().dropna(axis=0, how='all').dropna(axis=1, how='all')

    # This file type starts with metadata rows and then a 16-row x 24-column plate matrix.
    # Row labels like A, B, C, ... live in column 1 (or column 2 after removing leading NaN),
    # while the actual absorbance values are arranged in 3 blocks of 8 concentrations.
    for idx, row in df.iterrows():
        well = None
        cell_1 = row.iloc[0] if len(row) > 0 else None
        cell_2 = row.iloc[1] if len(row) > 1 else None

        if isinstance(cell_1, str) and len(cell_1.strip()) == 1 and cell_1.strip().isalpha():
            well = cell_1.strip()
        elif isinstance(cell_2, str) and len(cell_2.strip()) == 1 and cell_2.strip().isalpha():
            well = cell_2.strip()

        if well is None:
            continue

        # Extract the plate matrix values (columns after the row label)
        plate_values = []
        start_col = 1 if isinstance(cell_1, str) and len(cell_1.strip()) == 1 and cell_1.strip().isalpha() else 2
        for val in row.iloc[start_col: start_col + 24]:
            plate_values.append(parse_absorbance_value(val))

        if len(plate_values) < 24:
            continue

        for block_idx, ph_name in enumerate(PH_CONDITIONS.keys()):
            block_values = plate_values[block_idx * 8: (block_idx + 1) * 8]
            if len(block_values) != 8:
                continue
            for conc, abs_value in zip(DCPIP_CONCENTRATIONS, block_values):
                if not np.isfinite(abs_value):
                    continue
                data[ph_name]['conc'].append(conc)
                data[ph_name]['abs'].append(abs_value)

    # Fallback for legacy format with well labels in the first column (e.g. A1, A2)
    if all(len(v) == 0 for v in data.values() for v in [data[next(iter(data))]['conc'], data[next(iter(data))]['abs']]):
        print("Detected format: Well labels in first column")
        return parse_format_with_labels(df)

    print("Detected format: Matrix layout with row labels and three pH blocks")
    return data

def parse_format_with_labels(df):
    """Parse format where wells are labeled in the first column."""
    data = {ph: {'conc': [], 'abs': []} for ph in PH_CONDITIONS}
    
    for idx, row in df.iterrows():
        well = str(row.iloc[0]).strip()
        if len(well) >= 2 and well[0].isalpha():
            row_idx, col_idx = well_to_coordinates(well)
            col_num = col_idx + 1  # Convert to 1-indexed column number
            
            # Get absorbance value (usually in second column)
            abs_value = parse_absorbance_value(row.iloc[1] if len(row) > 1 else np.nan)
            if not np.isfinite(abs_value):
                continue
            
            # Determine which pH condition this well belongs to
            for ph_name, (col_start, col_end, _) in PH_CONDITIONS.items():
                if col_start <= col_num <= col_end:
                    # Map column to concentration (same pattern for each pH block)
                    conc_idx = (col_num - col_start) % 8
                    if conc_idx < len(DCPIP_CONCENTRATIONS):
                        conc = DCPIP_CONCENTRATIONS[conc_idx]
                        data[ph_name]['conc'].append(conc)
                        data[ph_name]['abs'].append(abs_value)
                    break
    
    return data

def parse_format_matrix(df):
    """Parse format where absorbance values are in a matrix (rows=wells, cols=replicates)."""
    data = {ph: {'conc': [], 'abs': []} for ph in PH_CONDITIONS}

    # The read-out for this file has row labels in a column and 24 numerical absorbance values
    # arranged as 3 pH blocks of 8 concentrations per row.
    for idx, row in df.iterrows():
        well = None
        if len(row) > 1 and isinstance(row.iloc[1], str) and len(row.iloc[1].strip()) == 1 and row.iloc[1].strip().isalpha():
            well = row.iloc[1].strip()
        elif len(row) > 0 and isinstance(row.iloc[0], str) and len(row.iloc[0].strip()) == 1 and row.iloc[0].strip().isalpha():
            well = row.iloc[0].strip()

        if well is None:
            continue

        values = []
        start_col = 2 if len(row) > 1 and isinstance(row.iloc[1], str) and len(row.iloc[1].strip()) == 1 and row.iloc[1].strip().isalpha() else 1
        for val in row.iloc[start_col:start_col + 24]:
            values.append(parse_absorbance_value(val))

        if len(values) < 24:
            continue

        for block_idx, ph_name in enumerate(PH_CONDITIONS.keys()):
            block_values = values[block_idx * 8: (block_idx + 1) * 8]
            if len(block_values) == 8:
                for conc, abs_value in zip(DCPIP_CONCENTRATIONS, block_values):
                    if not np.isfinite(abs_value):
                        continue
                    data[ph_name]['conc'].append(conc)
                    data[ph_name]['abs'].append(abs_value)

    return data

def organize_by_concentration(data):
    """
    Reorganize data by concentration to get mean absorbance for each concentration.
    """
    organized = {}
    for ph_name in PH_CONDITIONS:
        organized[ph_name] = {}
        conc_list = data[ph_name]['conc']
        abs_list = data[ph_name]['abs']
        
        for conc in set(conc_list):
            abs_values = [abs_list[i] for i in range(len(conc_list)) if conc_list[i] == conc and np.isfinite(abs_list[i])]
            if not abs_values:
                continue
            organized[ph_name][conc] = {
                'mean': np.mean(abs_values),
                'std': np.std(abs_values),
                'n': len(abs_values)
            }
    
    return organized

def plot_titration(organized_data, output_dir):
    """Create a single combined absorbance-vs-concentration plot with trendlines and summary metrics."""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle('DCPIP Titration Analysis', fontsize=16, fontweight='bold')

    for ph_name, (_, _, color) in PH_CONDITIONS.items():
        data_dict = organized_data[ph_name]

        conc_vals = sorted(data_dict.keys(), key=lambda x: x if x > 0 else float('inf'))
        abs_means = [data_dict[c]['mean'] for c in conc_vals]
        abs_stds = [data_dict[c]['std'] for c in conc_vals]

        if not conc_vals:
            continue

        x = np.asarray(conc_vals, dtype=float)
        y = np.asarray(abs_means, dtype=float)
        valid = np.isfinite(x) & np.isfinite(y) & (x >= 0)
        x_fit = x[valid]
        y_fit = y[valid]

        if x_fit.size >= 2:
            slope, intercept = np.polyfit(x_fit, y_fit, 1)
            y_pred = slope * x_fit + intercept
            ss_res = np.sum((y_fit - y_pred) ** 2)
            ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
            r2 = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot
            apparent_extinction = slope
            legend_label = (
                f'{ph_name}: ε_app = {apparent_extinction:.4g} AU·µM^-1, '
                f'R² = {r2:.3f}, slope = {slope:.4g}'
            )
            x_line = np.linspace(np.min(x_fit), np.max(x_fit), 200)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, linestyle='--', color=color, linewidth=1.8, alpha=0.8, zorder=2)
        else:
            legend_label = f'{ph_name}: insufficient data'

        ax.errorbar(
            x,
            y,
            yerr=abs_stds,
            fmt='o',
            linestyle='none',
            color=color,
            ecolor=color,
            elinewidth=1.2,
            capsize=3,
            markersize=6,
            alpha=0.95,
            zorder=3,
            label=legend_label
        )

    ax.set_xlabel('DCPIP Concentration (µM)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Absorbance (AU)', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim(0, 300)
    ax.set_ylim(bottom=0)
    ax.legend(loc='best', frameon=True, fontsize=9)

    plt.tight_layout()

    output_dir = Path(output_dir).expanduser()
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / 'DCPIP_titration_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nPlot saved to: {output_file}")
    except Exception as exc:
        fallback_dir = Path.cwd()
        output_file = fallback_dir / 'DCPIP_titration_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\nCould not save to {output_dir}. Saved instead to: {output_file}")
        print(f"Save error: {exc}")

    plt.show()
    return fig


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze DCPIP titration absorbance data from an Excel plate-reader export.'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        type=str,
        help='Path to the Excel file to analyze (for example: data.xlsx).'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Directory where the output plot will be saved. Defaults to the input file folder.'
    )
    return parser.parse_args()


def main():
    """Main analysis function."""
    args = parse_args()
    input_file = Path(args.input_file).expanduser() if args.input_file else None

    if input_file is None or not input_file.exists():
        print('ERROR: Please provide a valid Excel file path.')
        print('Example: python DCPIP_titrations_OT2_analysis.py "C:/path/to/file.xlsx"')
        print('Optional output directory: --output-dir "C:/path/to/output_folder"')
        return

    output_dir = Path(args.output_dir).expanduser() if args.output_dir else input_file.parent

    print("="*70)
    print("DCPIP TITRATION ANALYSIS")
    print("="*70)
    print(f"\nReading file: {input_file}")

    try:
        df = pd.read_excel(input_file)
        print("✓ File loaded successfully")
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns)}\n")
    except Exception as e:
        print(f"ERROR reading file: {e}")
        return

    print("Parsing plate data...")
    data = parse_plate_data(df)

    print("Organizing data by concentration...")
    organized = organize_by_concentration(data)

    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    for ph_name in PH_CONDITIONS:
        print(f"\n{ph_name}:")
        for conc in sorted(organized[ph_name].keys(), key=lambda x: x if x > 0 else float('inf')):
            info = organized[ph_name][conc]
            print(f"  {conc:7.2f} µM: Absorbance = {info['mean']:.4f} ± {info['std']:.4f} (n={int(info['n'])})")

    print("\n" + "="*70)
    print("Generating plots...")
    plot_titration(organized, output_dir)

    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()