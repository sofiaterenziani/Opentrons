"""
DCPIP Titrations Analysis
=========================
Reads a plate-reader Excel file from a 384-well plate DCPIP serial-dilution
experiment and produces a titration-curve figure.

Plate layout
------------
* Each reaction volume: 60 µL
* Titration step: 30 µL transferred into 60 µL → dilution factor = 30/90 = 1/3
  So each successive row has [DCPIP] = previous × (1/3)
* Starting concentration (row 1 of each section): 500 µM DCPIP
* Plate sections (by column in a 384-well plate, 24 columns total):
    Columns  1–8  : MES   pH 6
    Columns  9–16 : HEPES pH 7
    Columns 17–24 : HEPES pH 8

Usage
-----
    python DCPIP_titrations_analysis.py --file path/to/data.xlsx

The script expects the Excel file to contain absorbance values arranged so that
rows correspond to plate rows (A–P, i.e. 16 rows) and columns correspond to
plate columns (1–24).  If the sheet has a header row and/or row-label column
they are detected automatically.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import linregress


# ---------------------------------------------------------------------------
# Plate / dilution constants
# ---------------------------------------------------------------------------
PLATE_ROWS = 16          # A–P
PLATE_COLS = 24          # 1–24

INITIAL_DCPIP_uM = 500   # µM  (first row of each section)
TRANSFER_VOL_uL = 30     # µL transferred to next row
DEST_VOL_uL = 60         # µL already present in destination

DILUTION_FACTOR = TRANSFER_VOL_uL / (TRANSFER_VOL_uL + DEST_VOL_uL)  # = 1/3

SECTIONS = {
    "MES pH 6":   (0, 8),   # columns 0-7  (0-indexed)
    "HEPES pH 7": (8, 16),  # columns 8-15
    "HEPES pH 8": (16, 24), # columns 16-23
}

SECTION_COLORS = {
    "MES pH 6":   "#e6194b",
    "HEPES pH 7": "#3cb44b",
    "HEPES pH 8": "#4363d8",
}

SECTION_MARKERS = {
    "MES pH 6":   "o",
    "HEPES pH 7": "s",
    "HEPES pH 8": "^",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compute_dcpip_concentrations(n_steps: int) -> np.ndarray:
    """Return DCPIP concentrations (µM) for each titration step (row).

    Step 0 is the first row (undiluted, 500 µM).
    Each subsequent step is diluted by DILUTION_FACTOR.
    """
    steps = np.arange(n_steps)
    return INITIAL_DCPIP_uM * (DILUTION_FACTOR ** steps)


def load_plate_data(filepath: str) -> pd.DataFrame:
    """Load the Excel plate-reader file and return a (16 × 24) DataFrame.

    The function is tolerant of:
    * A leading row-label column (e.g. 'A', 'B', …)
    * A single header row with column numbers or names
    * Extra blank rows / columns
    """
    path = Path(filepath)
    if not path.exists():
        sys.exit(f"Error: file not found — {filepath}")

    raw = pd.read_excel(filepath, header=None)

    # Drop fully-empty rows and columns
    raw = raw.dropna(how="all", axis=0).dropna(how="all", axis=1).reset_index(drop=True)

    # Attempt to detect a header row: skip if the majority of values are non-numeric
    first_row = raw.iloc[0]
    numeric_mask = pd.to_numeric(first_row, errors="coerce").notna()
    if numeric_mask.sum() < len(first_row) / 2:
        raw = raw.iloc[1:].reset_index(drop=True)

    # Attempt to detect a row-label column: skip if the majority of values are non-numeric
    first_col = raw.iloc[:, 0]
    numeric_mask_col = pd.to_numeric(first_col, errors="coerce").notna()
    if numeric_mask_col.sum() < len(first_col) / 2:
        raw = raw.iloc[:, 1:].reset_index(drop=True)

    # Convert to numeric
    data = raw.apply(pd.to_numeric, errors="coerce")

    if data.shape[0] < PLATE_ROWS or data.shape[1] < PLATE_COLS:
        sys.exit(
            f"Error: expected at least {PLATE_ROWS} rows × {PLATE_COLS} columns "
            f"of numeric data, but found {data.shape[0]} × {data.shape[1]}."
        )

    data = data.iloc[:PLATE_ROWS, :PLATE_COLS]
    data.columns = range(PLATE_COLS)
    data.index = list("ABCDEFGHIJKLMNOP")
    return data


def section_stats(data: pd.DataFrame, col_start: int, col_end: int):
    """Return (mean, std) absorbance arrays per row for the given column slice."""
    section = data.iloc[:, col_start:col_end]
    return section.mean(axis=1).values, section.std(axis=1).values


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def make_figure(data: pd.DataFrame, output_path: str) -> None:
    concentrations = compute_dcpip_concentrations(PLATE_ROWS)  # µM, length 16

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    fig.suptitle("DCPIP Titration Curves — 384-Well Plate", fontsize=14, fontweight="bold")

    for ax, (section_name, (col_start, col_end)) in zip(axes, SECTIONS.items()):
        mean_abs, std_abs = section_stats(data, col_start, col_end)
        color = SECTION_COLORS[section_name]
        marker = SECTION_MARKERS[section_name]

        # Per-replicate (per-column) lines, semi-transparent
        for col_idx in range(col_start, col_end):
            ax.plot(concentrations, data.iloc[:, col_idx].values,
                    color=color, alpha=0.15, linewidth=0.8)

        # Mean ± std error bars
        ax.errorbar(concentrations, mean_abs, yerr=std_abs,
                    marker=marker, color=color, linewidth=1.8, markersize=6,
                    capsize=3, label="Mean ± SD")

        # Linear regression on the linear portion (first 8 points)
        n_linear = min(8, PLATE_ROWS)
        slope, intercept, r_value, p_value, _ = linregress(
            concentrations[:n_linear], mean_abs[:n_linear])
        x_fit = np.linspace(concentrations[n_linear - 1], concentrations[0], 200)
        y_fit = slope * x_fit + intercept
        # slope is AU/µM; multiply by 1e6 to convert to M⁻¹cm⁻¹ (Beer–Lambert, 1 cm path)
        epsilon = slope * 1e6
        ax.plot(x_fit, y_fit, "--", color="black", linewidth=1.2, alpha=0.7,
                label=f"Linear fit (R²={r_value**2:.3f})\nε={epsilon:.0f} M⁻¹cm⁻¹")

        ax.set_xscale("log")
        ax.set_xlabel("[DCPIP] (µM)", fontsize=11)
        ax.set_ylabel("Absorbance (a.u.)", fontsize=11)
        ax.set_title(section_name, fontsize=12)
        ax.legend(fontsize=8, loc="upper left")
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
        ax.grid(True, which="both", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Figure saved → {output_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Second figure: overlay of all three sections
# ---------------------------------------------------------------------------

def make_overlay_figure(data: pd.DataFrame, output_path: str) -> None:
    concentrations = compute_dcpip_concentrations(PLATE_ROWS)

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.suptitle("DCPIP Titration — All pH Conditions", fontsize=13, fontweight="bold")

    for section_name, (col_start, col_end) in SECTIONS.items():
        mean_abs, _ = section_stats(data, col_start, col_end)
        color = SECTION_COLORS[section_name]
        marker = SECTION_MARKERS[section_name]
        ax.plot(concentrations, mean_abs, marker=marker, color=color,
                linewidth=2, markersize=7, label=section_name)

    ax.set_xscale("log")
    ax.set_xlabel("[DCPIP] (µM)", fontsize=12)
    ax.set_ylabel("Absorbance (a.u.)", fontsize=12)
    ax.legend(fontsize=10)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:.1f}"))
    ax.grid(True, which="both", linestyle="--", alpha=0.4)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Overlay figure saved → {output_path}")
    plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyse DCPIP titration data from a 384-well plate reader Excel file."
    )
    parser.add_argument(
        "--file", "-f",
        required=True,
        help="Path to the plate-reader Excel file (.xlsx or .xls).",
    )
    parser.add_argument(
        "--out", "-o",
        default=None,
        help="Output figure file path (default: same directory as input, "
             "named DCPIP_titration_curves.png).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Determine output paths
    in_path = Path(args.file)
    if args.out:
        out_dir = Path(args.out).parent
        out_base = Path(args.out).stem
    else:
        out_dir = in_path.parent
        out_base = "DCPIP_titration_curves"

    out_path = str(out_dir / f"{out_base}.png")
    out_overlay = str(out_dir / f"{out_base}_overlay.png")

    print(f"Loading plate data from: {args.file}")
    data = load_plate_data(args.file)
    print(f"Data shape: {data.shape[0]} rows × {data.shape[1]} columns")

    # Print a summary of concentrations used
    concs = compute_dcpip_concentrations(PLATE_ROWS)
    print("\nDCPIP concentrations per row (µM):")
    row_labels = list("ABCDEFGHIJKLMNOP")
    for label, c in zip(row_labels, concs):
        print(f"  Row {label}: {c:.4f} µM")

    make_figure(data, out_path)
    make_overlay_figure(data, out_overlay)


if __name__ == "__main__":
    main()
