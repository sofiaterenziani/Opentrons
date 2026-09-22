"""Compare PedE and PedH kcat screens across pH, substrate, and metal descriptors.

Usage
-----
python compare_proteins.py \
    --condition PedE 6 path/to/PedE_pH6_avg_wells.csv \
    --condition PedE 7 path/to/PedE_pH7_avg_wells.csv \
    --condition PedE 8 path/to/PedE_pH8_avg_wells.csv \
    --condition PedH 6 path/to/PedH_pH6_avg_wells.csv \
    --condition PedH 7 path/to/PedH_pH7_avg_wells.csv \
    --condition PedH 8 path/to/PedH_pH8_avg_wells.csv \
    --outdir outputs/

The input CSVs may be single-plate files with ``kcat_s`` or replicate-averaged
files with ``kcat_mean_s`` / ``kcat_std_s``. The script writes long-form tables,
condition means, PedE-vs-PedH matched comparisons, descriptor correlation tables,
heatmaps, and descriptor scatter plots.

Notes
-----
1. Substrate and metal descriptors are bundled for the standard 15-substrate /
   24-metal screen and can be overridden with CSV files.
2. Predicted charge and logD are recalculated at each assay pH.
3. Apparent metal Kd or half-maximal metal concentration cannot be inferred from
   a single-concentration metal identity screen. A template CSV is exported so
   those values can be supplied later if separate titration data exist.
4. Oxidation potentials are left blank by default unless you provide them via the
   optional substrate descriptor override file.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import tempfile
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


METAL_ORDER = [
    "Al", "Ca", "Sc", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y",
    "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu",
]

CONTROL_SUBSTRATES = {"", "water", "tcep", "no_substrate"}
CONTROL_METALS = {"", "no_metal"}
DEFAULT_OUTDIR = "outputs"

NATURE_RC = {
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.title_fontsize": 7,
    "figure.titlesize": 7,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.0,
    "ytick.major.size": 2.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
}


@dataclass(frozen=True)
class SubstrateDescriptor:
    substrate: str
    substrate_class: str
    smiles: str
    carbon_number: int
    linear_series_name: str
    linear_series_index: int | None
    molecular_weight_da: float
    size_descriptor_name: str
    size_descriptor_value: float
    logp: float
    tpsa_a2: float
    rotatable_bonds: int
    hydroxyl_count: int
    alcohol_identity: str
    aromatic_ring_count: int
    acid_pka: float | None = None
    base_pka: float | None = None
    oxidation_potential_v: float | None = None
    oxidation_potential_source: str = ""
    notes: str = ""


@dataclass(frozen=True)
class MetalDescriptor:
    metal: str
    atomic_number: int
    charge: int
    ionic_radius_pm: float
    coordination_number: int
    f_electron_count: int
    is_lanthanide: bool
    source: str
    notes: str = ""


@dataclass
class WellRecord:
    protein: str
    pH: float
    well: str
    row: str
    col: int
    metal: str
    substrate: str
    kcat: float
    kcat_std: float


@dataclass
class ConditionMean:
    protein: str
    pH: float
    metal: str
    substrate: str
    mean_kcat: float
    std_kcat: float
    sem_kcat: float
    n_values: int
    metal_desc: MetalDescriptor
    substrate_desc: SubstrateDescriptor
    predicted_charge: float
    logd: float
    apparent_halfmax_metal_uM: float | None
    apparent_metal_kd_uM: float | None
    metal_kinetics_source: str


@dataclass(frozen=True)
class MetalKinetics:
    protein: str
    pH: float
    metal: str
    apparent_halfmax_metal_uM: float | None = None
    apparent_metal_kd_uM: float | None = None
    source: str = ""


DEFAULT_SUBSTRATES = {
    "methanol": SubstrateDescriptor(
        substrate="methanol", substrate_class="alcohol", smiles="CO", carbon_number=1,
        linear_series_name="linear_primary_alcohols", linear_series_index=1,
        molecular_weight_da=32.04, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=32.04, logp=-0.74, tpsa_a2=20.23, rotatable_bonds=0,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "ethanol": SubstrateDescriptor(
        substrate="ethanol", substrate_class="alcohol", smiles="CCO", carbon_number=2,
        linear_series_name="linear_primary_alcohols", linear_series_index=2,
        molecular_weight_da=46.07, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=46.07, logp=-0.31, tpsa_a2=20.23, rotatable_bonds=0,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "1-propanol": SubstrateDescriptor(
        substrate="1-propanol", substrate_class="alcohol", smiles="CCCO", carbon_number=3,
        linear_series_name="linear_primary_alcohols", linear_series_index=3,
        molecular_weight_da=60.10, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=60.10, logp=0.25, tpsa_a2=20.23, rotatable_bonds=1,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "1-butanol": SubstrateDescriptor(
        substrate="1-butanol", substrate_class="alcohol", smiles="CCCCO", carbon_number=4,
        linear_series_name="linear_primary_alcohols", linear_series_index=4,
        molecular_weight_da=74.12, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=74.12, logp=0.88, tpsa_a2=20.23, rotatable_bonds=2,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "1-pentanol": SubstrateDescriptor(
        substrate="1-pentanol", substrate_class="alcohol", smiles="CCCCCO", carbon_number=5,
        linear_series_name="linear_primary_alcohols", linear_series_index=5,
        molecular_weight_da=88.15, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=88.15, logp=1.51, tpsa_a2=20.23, rotatable_bonds=3,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "1-hexanol": SubstrateDescriptor(
        substrate="1-hexanol", substrate_class="alcohol", smiles="CCCCCCO", carbon_number=6,
        linear_series_name="linear_primary_alcohols", linear_series_index=6,
        molecular_weight_da=102.17, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=102.17, logp=2.03, tpsa_a2=20.23, rotatable_bonds=4,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "3-methyl-1-butanol": SubstrateDescriptor(
        substrate="3-methyl-1-butanol", substrate_class="alcohol", smiles="CC(C)CCO", carbon_number=5,
        linear_series_name="branched_primary_alcohols", linear_series_index=1,
        molecular_weight_da=88.15, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=88.15, logp=1.16, tpsa_a2=20.23, rotatable_bonds=2,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "2-methyl-1-butanol": SubstrateDescriptor(
        substrate="2-methyl-1-butanol", substrate_class="alcohol", smiles="CCC(C)CO", carbon_number=5,
        linear_series_name="branched_primary_alcohols", linear_series_index=2,
        molecular_weight_da=88.15, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=88.15, logp=1.16, tpsa_a2=20.23, rotatable_bonds=2,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "1,5-pentanediol": SubstrateDescriptor(
        substrate="1,5-pentanediol", substrate_class="alcohol", smiles="OCCCCCO", carbon_number=5,
        linear_series_name="diols", linear_series_index=1,
        molecular_weight_da=104.15, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=104.15, logp=0.14, tpsa_a2=40.46, rotatable_bonds=4,
        hydroxyl_count=2, alcohol_identity="primary_diol", aromatic_ring_count=0,
        notes="Built-in fallback descriptors"
    ),
    "2-phenylethanol": SubstrateDescriptor(
        substrate="2-phenylethanol", substrate_class="alcohol", smiles="OCCc1ccccc1", carbon_number=8,
        linear_series_name="aryl_alcohols", linear_series_index=1,
        molecular_weight_da=122.16, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=122.16, logp=1.36, tpsa_a2=20.23, rotatable_bonds=2,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=1,
        notes="Built-in fallback descriptors"
    ),
    "4-hydroxybenzyl alcohol": SubstrateDescriptor(
        substrate="4-hydroxybenzyl alcohol", substrate_class="alcohol", smiles="OCC1=CC=C(O)C=C1", carbon_number=7,
        linear_series_name="aryl_alcohols", linear_series_index=2,
        molecular_weight_da=124.14, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=124.14, logp=0.73, tpsa_a2=40.46, rotatable_bonds=1,
        hydroxyl_count=2, alcohol_identity="primary", aromatic_ring_count=1,
        notes="Built-in fallback descriptors"
    ),
    "vanillin": SubstrateDescriptor(
        substrate="vanillin", substrate_class="aldehyde", smiles="COc1cc(C=O)ccc1O", carbon_number=8,
        linear_series_name="aryl_aldehydes", linear_series_index=1,
        molecular_weight_da=152.15, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=152.15, logp=1.21, tpsa_a2=46.53, rotatable_bonds=2,
        hydroxyl_count=1, alcohol_identity="not_alcohol", aromatic_ring_count=1,
        notes="Built-in fallback descriptors"
    ),
    "vanillyl alcohol": SubstrateDescriptor(
        substrate="vanillyl alcohol", substrate_class="alcohol", smiles="COc1cc(CO)ccc1O", carbon_number=8,
        linear_series_name="aryl_alcohols", linear_series_index=3,
        molecular_weight_da=154.16, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=154.16, logp=0.85, tpsa_a2=49.69, rotatable_bonds=2,
        hydroxyl_count=2, alcohol_identity="primary", aromatic_ring_count=1,
        notes="Built-in fallback descriptors"
    ),
    "protocatechuic acid": SubstrateDescriptor(
        substrate="protocatechuic acid", substrate_class="acid", smiles="O=C(O)c1ccc(O)c(O)c1", carbon_number=7,
        linear_series_name="aryl_acids", linear_series_index=1,
        molecular_weight_da=154.12, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=154.12, logp=1.10, tpsa_a2=77.76, rotatable_bonds=1,
        hydroxyl_count=2, alcohol_identity="not_alcohol", aromatic_ring_count=1,
        acid_pka=4.00, notes="Built-in fallback descriptors"
    ),
    "vanillic acid": SubstrateDescriptor(
        substrate="vanillic acid", substrate_class="acid", smiles="COc1cc(C(=O)O)ccc1O", carbon_number=8,
        linear_series_name="aryl_acids", linear_series_index=2,
        molecular_weight_da=168.15, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=168.15, logp=1.24, tpsa_a2=66.76, rotatable_bonds=2,
        hydroxyl_count=1, alcohol_identity="not_alcohol", aromatic_ring_count=1,
        acid_pka=4.16, notes="Built-in fallback descriptors"
    ),
    "ethanolamine": SubstrateDescriptor(
        substrate="ethanolamine", substrate_class="amino_alcohol", smiles="NCCO", carbon_number=2,
        linear_series_name="amino_alcohols", linear_series_index=1,
        molecular_weight_da=61.08, size_descriptor_name="molecular_weight_da",
        size_descriptor_value=61.08, logp=-1.31, tpsa_a2=46.25, rotatable_bonds=1,
        hydroxyl_count=1, alcohol_identity="primary", aromatic_ring_count=0,
        base_pka=9.50, notes="Built-in fallback descriptors"
    ),
}


DEFAULT_METALS = {
    "Al": MetalDescriptor("Al", 13, 3, 53.5, 6, 0, False, "Shannon CN=6 radius"),
    "Ca": MetalDescriptor("Ca", 20, 2, 100.0, 6, 0, False, "Shannon CN=6 radius"),
    "Sc": MetalDescriptor("Sc", 21, 3, 74.5, 6, 0, False, "Shannon CN=6 radius"),
    "Mn": MetalDescriptor("Mn", 25, 2, 83.0, 6, 0, False, "Shannon CN=6 radius; Mn2+"),
    "Fe": MetalDescriptor("Fe", 26, 3, 64.5, 6, 0, False, "Shannon CN=6 radius; Fe3+"),
    "Co": MetalDescriptor("Co", 27, 2, 74.5, 6, 0, False, "Shannon CN=6 radius; Co2+"),
    "Ni": MetalDescriptor("Ni", 28, 2, 69.0, 6, 0, False, "Shannon CN=6 radius; Ni2+"),
    "Cu": MetalDescriptor("Cu", 29, 2, 73.0, 6, 0, False, "Shannon CN=6 radius; Cu2+"),
    "Zn": MetalDescriptor("Zn", 30, 2, 74.0, 6, 0, False, "Shannon CN=6 radius; Zn2+"),
    "Y": MetalDescriptor("Y", 39, 3, 90.0, 6, 0, False, "Shannon CN=6 radius"),
    "La": MetalDescriptor("La", 57, 3, 103.2, 6, 0, True, "Shannon CN=6 radius"),
    "Ce": MetalDescriptor("Ce", 58, 3, 101.0, 6, 1, True, "Shannon CN=6 radius"),
    "Pr": MetalDescriptor("Pr", 59, 3, 99.0, 6, 2, True, "Shannon CN=6 radius"),
    "Nd": MetalDescriptor("Nd", 60, 3, 98.3, 6, 3, True, "Shannon CN=6 radius"),
    "Sm": MetalDescriptor("Sm", 62, 3, 95.8, 6, 5, True, "Shannon CN=6 radius"),
    "Eu": MetalDescriptor("Eu", 63, 3, 94.7, 6, 6, True, "Shannon CN=6 radius"),
    "Gd": MetalDescriptor("Gd", 64, 3, 93.8, 6, 7, True, "Shannon CN=6 radius"),
    "Tb": MetalDescriptor("Tb", 65, 3, 92.3, 6, 8, True, "Shannon CN=6 radius"),
    "Dy": MetalDescriptor("Dy", 66, 3, 91.2, 6, 9, True, "Shannon CN=6 radius"),
    "Ho": MetalDescriptor("Ho", 67, 3, 90.1, 6, 10, True, "Shannon CN=6 radius"),
    "Er": MetalDescriptor("Er", 68, 3, 89.0, 6, 11, True, "Shannon CN=6 radius"),
    "Tm": MetalDescriptor("Tm", 69, 3, 88.0, 6, 12, True, "Shannon CN=6 radius"),
    "Yb": MetalDescriptor("Yb", 70, 3, 86.8, 6, 13, True, "Shannon CN=6 radius"),
    "Lu": MetalDescriptor("Lu", 71, 3, 86.1, 6, 14, True, "Shannon CN=6 radius"),
}


def normalize_name(value: str) -> str:
    return " ".join(value.strip().lower().split())


def finite_or_blank(value: float | int | str | None) -> str | float | int:
    if value is None:
        return ""
    if isinstance(value, float) and not np.isfinite(value):
        return ""
    return value


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--manifest", type=Path, default=None,
                        help="JSON manifest with {conditions: [{protein, pH, csv}, ...]}.")
    parser.add_argument("--condition", nargs=3, action="append", metavar=("PROTEIN", "PH", "CSV"),
                        help="Condition triplet: PROTEIN PH CSV. May be repeated.")
    parser.add_argument("--outdir", type=Path, default=Path(DEFAULT_OUTDIR),
                        help="Output directory for CSVs and figures.")
    parser.add_argument("--label", type=str, default="compare",
                        help="Prefix for output file names.")
    parser.add_argument("--substrate-descriptors", type=Path, default=None,
                        help="Optional CSV overriding built-in substrate descriptors.")
    parser.add_argument("--metal-descriptors", type=Path, default=None,
                        help="Optional CSV overriding built-in metal descriptors.")
    parser.add_argument("--metal-kinetics", type=Path, default=None,
                        help="Optional CSV with apparent half-max metal concentration or Kd values.")
    parser.add_argument("--kcat-upper", type=float, default=None,
                        help="Optional upper y-limit / color limit for plots.")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip figure generation and write CSVs only.")
    return parser


def resolve_conditions(args: argparse.Namespace, parser: argparse.ArgumentParser) -> list[tuple[str, float, Path]]:
    if args.manifest is not None:
        payload = json.loads(args.manifest.read_text())
        raw_items = payload.get("conditions", payload)
        conditions = []
        for item in raw_items:
            try:
                protein = str(item["protein"])
                pH = float(item["pH"])
                csv_path = Path(item["csv"])
            except (TypeError, ValueError, KeyError) as exc:
                parser.error(f"Malformed manifest entry {item!r}: {exc}")
            if not csv_path.is_absolute():
                csv_path = (args.manifest.parent / csv_path).resolve()
            conditions.append((protein, pH, csv_path))
        if not conditions:
            parser.error("Manifest did not contain any conditions.")
        return conditions
    if not args.condition:
        parser.error("Provide either --manifest or one or more --condition PROTEIN PH CSV entries.")
    conditions = []
    for protein, ph_text, csv_text in args.condition:
        try:
            pH = float(ph_text)
        except ValueError:
            parser.error(f"Invalid pH value {ph_text!r}.")
        conditions.append((protein, pH, Path(csv_text).resolve()))
    return conditions


def read_wells_csv(csv_path: Path, protein: str, pH: float) -> list[WellRecord]:
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        return []
    first = rows[0]
    if "kcat_mean_s" in first:
        kcat_key = "kcat_mean_s"
        std_key = "kcat_std_s" if "kcat_std_s" in first else None
    elif "kcat_s" in first:
        kcat_key = "kcat_s"
        std_key = None
    else:
        raise ValueError(f"{csv_path}: expected kcat_mean_s or kcat_s column.")
    records = []
    for row in rows:
        substrate = row.get("substrate", "").strip()
        metal = row.get("metal", "").strip()
        if normalize_name(substrate) in CONTROL_SUBSTRATES or normalize_name(metal) in CONTROL_METALS:
            continue
        raw_kcat = row.get(kcat_key, "")
        if raw_kcat in ("", "nan", "NaN"):
            continue
        try:
            kcat = float(raw_kcat)
        except ValueError:
            continue
        if not np.isfinite(kcat):
            continue
        kcat_std = float("nan")
        if std_key is not None and row.get(std_key, "") not in ("", "nan", "NaN"):
            try:
                kcat_std = float(row[std_key])
            except ValueError:
                kcat_std = float("nan")
        records.append(WellRecord(
            protein=protein,
            pH=pH,
            well=row.get("well", ""),
            row=row.get("row", ""),
            col=int(row.get("col", 0) or 0),
            metal=metal,
            substrate=substrate,
            kcat=kcat,
            kcat_std=kcat_std,
        ))
    return records


def parse_optional_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if text in ("", "nan", "NaN", "none", "None"):
        return None
    return float(text)


def load_substrate_descriptors(path: Path | None) -> dict[str, SubstrateDescriptor]:
    descriptors = dict(DEFAULT_SUBSTRATES)
    if path is None:
        return descriptors
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {
            "substrate", "substrate_class", "smiles", "carbon_number", "linear_series_name",
            "linear_series_index", "molecular_weight_da", "size_descriptor_name",
            "size_descriptor_value", "logp", "tpsa_a2", "rotatable_bonds",
            "hydroxyl_count", "alcohol_identity", "aromatic_ring_count",
            "acid_pka", "base_pka", "oxidation_potential_v", "oxidation_potential_source", "notes",
        }
        missing_columns = required - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(f"{path}: missing columns {sorted(missing_columns)}")
        for row in reader:
            descriptor = SubstrateDescriptor(
                substrate=row["substrate"].strip(),
                substrate_class=row["substrate_class"].strip(),
                smiles=row["smiles"].strip(),
                carbon_number=int(row["carbon_number"]),
                linear_series_name=row["linear_series_name"].strip(),
                linear_series_index=int(row["linear_series_index"]) if row["linear_series_index"].strip() else None,
                molecular_weight_da=float(row["molecular_weight_da"]),
                size_descriptor_name=row["size_descriptor_name"].strip(),
                size_descriptor_value=float(row["size_descriptor_value"]),
                logp=float(row["logp"]),
                tpsa_a2=float(row["tpsa_a2"]),
                rotatable_bonds=int(row["rotatable_bonds"]),
                hydroxyl_count=int(row["hydroxyl_count"]),
                alcohol_identity=row["alcohol_identity"].strip(),
                aromatic_ring_count=int(row["aromatic_ring_count"]),
                acid_pka=parse_optional_float(row.get("acid_pka")),
                base_pka=parse_optional_float(row.get("base_pka")),
                oxidation_potential_v=parse_optional_float(row.get("oxidation_potential_v")),
                oxidation_potential_source=row.get("oxidation_potential_source", "").strip(),
                notes=row.get("notes", "").strip(),
            )
            descriptors[normalize_name(descriptor.substrate)] = descriptor
    return descriptors


def load_metal_descriptors(path: Path | None) -> dict[str, MetalDescriptor]:
    descriptors = dict(DEFAULT_METALS)
    if path is None:
        return descriptors
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {
            "metal", "atomic_number", "charge", "ionic_radius_pm", "coordination_number",
            "f_electron_count", "is_lanthanide", "source", "notes",
        }
        missing_columns = required - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(f"{path}: missing columns {sorted(missing_columns)}")
        for row in reader:
            descriptor = MetalDescriptor(
                metal=row["metal"].strip(),
                atomic_number=int(row["atomic_number"]),
                charge=int(row["charge"]),
                ionic_radius_pm=float(row["ionic_radius_pm"]),
                coordination_number=int(row["coordination_number"]),
                f_electron_count=int(row["f_electron_count"]),
                is_lanthanide=row["is_lanthanide"].strip().lower() in {"1", "true", "yes"},
                source=row["source"].strip(),
                notes=row.get("notes", "").strip(),
            )
            descriptors[descriptor.metal] = descriptor
    return descriptors


def load_metal_kinetics(path: Path | None) -> dict[tuple[str, float, str], MetalKinetics]:
    result: dict[tuple[str, float, str], MetalKinetics] = {}
    if path is None:
        return result
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"protein", "pH", "metal", "apparent_halfmax_metal_uM", "apparent_metal_kd_uM", "source"}
        missing_columns = required - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(f"{path}: missing columns {sorted(missing_columns)}")
        for row in reader:
            kinetics = MetalKinetics(
                protein=row["protein"].strip(),
                pH=float(row["pH"]),
                metal=row["metal"].strip(),
                apparent_halfmax_metal_uM=parse_optional_float(row.get("apparent_halfmax_metal_uM")),
                apparent_metal_kd_uM=parse_optional_float(row.get("apparent_metal_kd_uM")),
                source=row.get("source", "").strip(),
            )
            result[(kinetics.protein, kinetics.pH, kinetics.metal)] = kinetics
    return result


def substrate_charge(descriptor: SubstrateDescriptor, pH: float) -> float:
    charge = 0.0
    if descriptor.base_pka is not None:
        protonated_fraction = 1.0 / (1.0 + 10.0 ** (pH - descriptor.base_pka))
        charge += protonated_fraction
    if descriptor.acid_pka is not None:
        deprotonated_fraction = 1.0 / (1.0 + 10.0 ** (descriptor.acid_pka - pH))
        charge -= deprotonated_fraction
    return charge


def substrate_logd(descriptor: SubstrateDescriptor, pH: float) -> float:
    neutral_fraction = 1.0
    if descriptor.base_pka is not None:
        neutral_fraction *= 1.0 / (1.0 + 10.0 ** (descriptor.base_pka - pH))
    if descriptor.acid_pka is not None:
        neutral_fraction *= 1.0 / (1.0 + 10.0 ** (pH - descriptor.acid_pka))
    neutral_fraction = max(neutral_fraction, 1e-12)
    return descriptor.logp + math.log10(neutral_fraction)


def group_stats(values: list[float]) -> tuple[float, float, float, int]:
    array = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if array.size == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(np.mean(array))
    std = float(np.std(array, ddof=1)) if array.size > 1 else 0.0
    sem = float(std / np.sqrt(array.size)) if array.size > 1 else 0.0
    return mean, std, sem, int(array.size)


def resolve_substrate_descriptor(name: str, descriptors: dict[str, SubstrateDescriptor]) -> SubstrateDescriptor | None:
    key = normalize_name(name)
    return descriptors.get(key)


def aggregate_condition_means(
    records: list[WellRecord],
    substrate_descriptors: dict[str, SubstrateDescriptor],
    metal_descriptors: dict[str, MetalDescriptor],
    metal_kinetics: dict[tuple[str, float, str], MetalKinetics],
) -> list[ConditionMean]:
    grouped: dict[tuple[str, float, str, str], list[WellRecord]] = defaultdict(list)
    for record in records:
        grouped[(record.protein, record.pH, record.metal, record.substrate)].append(record)
    rows: list[ConditionMean] = []
    missing_substrates: set[str] = set()
    missing_metals: set[str] = set()
    for key in sorted(grouped):
        protein, pH, metal, substrate = key
        substrate_desc = resolve_substrate_descriptor(substrate, substrate_descriptors)
        metal_desc = metal_descriptors.get(metal)
        if substrate_desc is None:
            missing_substrates.add(substrate)
            continue
        if metal_desc is None:
            missing_metals.add(metal)
            continue
        mean, std, sem, n_values = group_stats([item.kcat for item in grouped[key]])
        kinetics = metal_kinetics.get((protein, pH, metal), MetalKinetics(protein, pH, metal))
        rows.append(ConditionMean(
            protein=protein,
            pH=pH,
            metal=metal,
            substrate=substrate,
            mean_kcat=mean,
            std_kcat=std,
            sem_kcat=sem,
            n_values=n_values,
            metal_desc=metal_desc,
            substrate_desc=replace(substrate_desc, substrate=substrate),
            predicted_charge=substrate_charge(substrate_desc, pH),
            logd=substrate_logd(substrate_desc, pH),
            apparent_halfmax_metal_uM=kinetics.apparent_halfmax_metal_uM,
            apparent_metal_kd_uM=kinetics.apparent_metal_kd_uM,
            metal_kinetics_source=kinetics.source,
        ))
    if missing_substrates or missing_metals:
        issues = []
        if missing_substrates:
            issues.append(f"missing substrate descriptors: {sorted(missing_substrates)}")
        if missing_metals:
            issues.append(f"missing metal descriptors: {sorted(missing_metals)}")
        raise ValueError("; ".join(issues))
    return rows


def write_wells_long(
    out_csv: Path,
    records: list[WellRecord],
    substrate_descriptors: dict[str, SubstrateDescriptor],
    metal_descriptors: dict[str, MetalDescriptor],
) -> None:
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "protein", "pH", "well", "row", "col", "metal", "substrate", "kcat_s", "kcat_std_s",
            "substrate_class", "carbon_number", "molecular_weight_da", "size_descriptor_name",
            "size_descriptor_value", "logp", "predicted_charge", "logd", "tpsa_a2",
            "rotatable_bonds", "hydroxyl_count", "alcohol_identity", "aromatic_ring_count",
            "acid_pka", "base_pka", "oxidation_potential_v", "metal_charge", "ionic_radius_pm",
            "coordination_number", "atomic_number", "f_electron_count", "charge_density_charge_per_a2",
        ])
        for record in records:
            substrate_desc = resolve_substrate_descriptor(record.substrate, substrate_descriptors)
            metal_desc = metal_descriptors.get(record.metal)
            if substrate_desc is None or metal_desc is None:
                continue
            writer.writerow([
                record.protein, record.pH, record.well, record.row, record.col, record.metal, record.substrate,
                finite_or_blank(record.kcat), finite_or_blank(record.kcat_std), substrate_desc.substrate_class,
                substrate_desc.carbon_number, substrate_desc.molecular_weight_da, substrate_desc.size_descriptor_name,
                substrate_desc.size_descriptor_value, substrate_desc.logp,
                substrate_charge(substrate_desc, record.pH), substrate_logd(substrate_desc, record.pH),
                substrate_desc.tpsa_a2, substrate_desc.rotatable_bonds, substrate_desc.hydroxyl_count,
                substrate_desc.alcohol_identity, substrate_desc.aromatic_ring_count,
                finite_or_blank(substrate_desc.acid_pka), finite_or_blank(substrate_desc.base_pka),
                finite_or_blank(substrate_desc.oxidation_potential_v), metal_desc.charge, metal_desc.ionic_radius_pm,
                metal_desc.coordination_number, metal_desc.atomic_number, metal_desc.f_electron_count,
                charge_density(metal_desc),
            ])


def write_condition_means(out_csv: Path, rows: list[ConditionMean]) -> None:
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "protein", "pH", "metal", "substrate", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_values",
            "substrate_class", "carbon_number", "molecular_weight_da", "size_descriptor_name",
            "size_descriptor_value", "logp", "predicted_charge", "logd", "tpsa_a2", "rotatable_bonds",
            "hydroxyl_count", "alcohol_identity", "aromatic_ring_count", "acid_pka", "base_pka",
            "oxidation_potential_v", "oxidation_potential_source", "metal_charge", "ionic_radius_pm",
            "coordination_number", "atomic_number", "f_electron_count", "charge_density_charge_per_a2",
            "is_lanthanide", "apparent_halfmax_metal_uM", "apparent_metal_kd_uM", "metal_kinetics_source",
        ])
        for row in rows:
            writer.writerow([
                row.protein, row.pH, row.metal, row.substrate, row.mean_kcat, row.std_kcat, row.sem_kcat, row.n_values,
                row.substrate_desc.substrate_class, row.substrate_desc.carbon_number,
                row.substrate_desc.molecular_weight_da, row.substrate_desc.size_descriptor_name,
                row.substrate_desc.size_descriptor_value, row.substrate_desc.logp, row.predicted_charge, row.logd,
                row.substrate_desc.tpsa_a2, row.substrate_desc.rotatable_bonds, row.substrate_desc.hydroxyl_count,
                row.substrate_desc.alcohol_identity, row.substrate_desc.aromatic_ring_count,
                finite_or_blank(row.substrate_desc.acid_pka), finite_or_blank(row.substrate_desc.base_pka),
                finite_or_blank(row.substrate_desc.oxidation_potential_v), row.substrate_desc.oxidation_potential_source,
                row.metal_desc.charge, row.metal_desc.ionic_radius_pm, row.metal_desc.coordination_number,
                row.metal_desc.atomic_number, row.metal_desc.f_electron_count, charge_density(row.metal_desc),
                int(row.metal_desc.is_lanthanide), finite_or_blank(row.apparent_halfmax_metal_uM),
                finite_or_blank(row.apparent_metal_kd_uM), row.metal_kinetics_source,
            ])


def charge_density(metal_desc: MetalDescriptor) -> float:
    radius_a = metal_desc.ionic_radius_pm / 100.0
    return metal_desc.charge / (radius_a ** 2)


def aggregate_by_substrate(rows: list[ConditionMean]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float, str], list[ConditionMean]] = defaultdict(list)
    for row in rows:
        grouped[(row.protein, row.pH, row.substrate)].append(row)
    summaries = []
    for (protein, pH, substrate), items in sorted(grouped.items()):
        mean, std, sem, n_values = group_stats([item.mean_kcat for item in items])
        exemplar = items[0]
        summaries.append({
            "protein": protein,
            "pH": pH,
            "substrate": substrate,
            "mean_kcat_s": mean,
            "std_kcat_s": std,
            "sem_kcat_s": sem,
            "n_metals": n_values,
            "substrate_class": exemplar.substrate_desc.substrate_class,
            "carbon_number": exemplar.substrate_desc.carbon_number,
            "molecular_weight_da": exemplar.substrate_desc.molecular_weight_da,
            "size_descriptor_name": exemplar.substrate_desc.size_descriptor_name,
            "size_descriptor_value": exemplar.substrate_desc.size_descriptor_value,
            "logp": exemplar.substrate_desc.logp,
            "predicted_charge": exemplar.predicted_charge,
            "logd": exemplar.logd,
            "tpsa_a2": exemplar.substrate_desc.tpsa_a2,
            "rotatable_bonds": exemplar.substrate_desc.rotatable_bonds,
            "hydroxyl_count": exemplar.substrate_desc.hydroxyl_count,
            "alcohol_identity": exemplar.substrate_desc.alcohol_identity,
            "aromatic_ring_count": exemplar.substrate_desc.aromatic_ring_count,
            "oxidation_potential_v": exemplar.substrate_desc.oxidation_potential_v,
            "linear_series_name": exemplar.substrate_desc.linear_series_name,
            "linear_series_index": exemplar.substrate_desc.linear_series_index,
        })
    return summaries


def aggregate_by_metal(rows: list[ConditionMean]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float, str], list[ConditionMean]] = defaultdict(list)
    for row in rows:
        grouped[(row.protein, row.pH, row.metal)].append(row)
    summaries = []
    for (protein, pH, metal), items in sorted(grouped.items()):
        mean, std, sem, n_values = group_stats([item.mean_kcat for item in items])
        exemplar = items[0]
        summaries.append({
            "protein": protein,
            "pH": pH,
            "metal": metal,
            "mean_kcat_s": mean,
            "std_kcat_s": std,
            "sem_kcat_s": sem,
            "n_substrates": n_values,
            "metal_charge": exemplar.metal_desc.charge,
            "ionic_radius_pm": exemplar.metal_desc.ionic_radius_pm,
            "coordination_number": exemplar.metal_desc.coordination_number,
            "atomic_number": exemplar.metal_desc.atomic_number,
            "f_electron_count": exemplar.metal_desc.f_electron_count,
            "charge_density_charge_per_a2": charge_density(exemplar.metal_desc),
            "is_lanthanide": exemplar.metal_desc.is_lanthanide,
            "apparent_halfmax_metal_uM": exemplar.apparent_halfmax_metal_uM,
            "apparent_metal_kd_uM": exemplar.apparent_metal_kd_uM,
            "metal_kinetics_source": exemplar.metal_kinetics_source,
        })
    return summaries


def write_dict_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            clean = {key: finite_or_blank(row.get(key)) for key in fieldnames}
            writer.writerow(clean)


def matched_protein_rows(rows: list[ConditionMean]) -> list[dict[str, object]]:
    proteins = sorted({row.protein for row in rows})
    if len(proteins) != 2:
        return []
    left, right = proteins
    lookup = {(row.protein, row.pH, row.metal, row.substrate): row for row in rows}
    matches = []
    for pH, metal, substrate in sorted({(row.pH, row.metal, row.substrate) for row in rows}):
        left_row = lookup.get((left, pH, metal, substrate))
        right_row = lookup.get((right, pH, metal, substrate))
        if left_row is None or right_row is None:
            continue
        log2_ratio = float("nan")
        if left_row.mean_kcat > 0 and right_row.mean_kcat > 0:
            log2_ratio = float(np.log2(right_row.mean_kcat / left_row.mean_kcat))
        matches.append({
            "pH": pH,
            "metal": metal,
            "substrate": substrate,
            f"{left}_kcat_s": left_row.mean_kcat,
            f"{right}_kcat_s": right_row.mean_kcat,
            f"{right}_minus_{left}_kcat_s": right_row.mean_kcat - left_row.mean_kcat,
            f"{right}_over_{left}_log2_ratio": log2_ratio,
            "substrate_class": left_row.substrate_desc.substrate_class,
            "carbon_number": left_row.substrate_desc.carbon_number,
            "logd": left_row.logd,
            "ionic_radius_pm": left_row.metal_desc.ionic_radius_pm,
            "f_electron_count": left_row.metal_desc.f_electron_count,
        })
    return matches


def finite_pairs(xs: list[object], ys: list[object]) -> tuple[np.ndarray, np.ndarray]:
    x_values = []
    y_values = []
    for x, y in zip(xs, ys):
        if x is None or y is None:
            continue
        try:
            x_float = float(x)
            y_float = float(y)
        except (TypeError, ValueError):
            continue
        if np.isfinite(x_float) and np.isfinite(y_float):
            x_values.append(x_float)
            y_values.append(y_float)
    return np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float)


def safe_correlation(xs: list[object], ys: list[object]) -> tuple[float, float, float, int]:
    x_array, y_array = finite_pairs(xs, ys)
    if x_array.size < 2:
        return float("nan"), float("nan"), float("nan"), int(x_array.size)
    if np.allclose(x_array, x_array[0]) or np.allclose(y_array, y_array[0]):
        return float("nan"), float("nan"), float("nan"), int(x_array.size)
    slope, intercept = np.polyfit(x_array, y_array, 1)
    corr = float(np.corrcoef(x_array, y_array)[0, 1])
    return corr, float(slope), float(intercept), int(x_array.size)


def descriptor_correlations(rows: list[dict[str, object]], key_name: str, descriptor_names: list[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["protein"]), float(row["pH"]))].append(row)
    results = []
    for (protein, pH), items in sorted(grouped.items()):
        y_values = [item["mean_kcat_s"] for item in items]
        for descriptor in descriptor_names:
            corr, slope, intercept, n_points = safe_correlation([item[descriptor] for item in items], y_values)
            results.append({
                "protein": protein,
                "pH": pH,
                "axis": key_name,
                "descriptor": descriptor,
                "pearson_r": corr,
                "slope": slope,
                "intercept": intercept,
                "n_points": n_points,
            })
    return results


def subset_correlations(rows: list[dict[str, object]], subset_name: str, subset_filter, descriptor_names: list[str]) -> list[dict[str, object]]:
    filtered = [row for row in rows if subset_filter(row)]
    if not filtered:
        return []
    results = descriptor_correlations(filtered, subset_name, descriptor_names)
    for row in results:
        row["subset"] = subset_name
    return results


def categorical_group_summary(rows: list[dict[str, object]], category_key: str) -> list[dict[str, object]]:
    grouped: dict[tuple[str, float, str], list[float]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["protein"]), float(row["pH"]), str(row[category_key]))].append(float(row["mean_kcat_s"]))
    summaries = []
    for (protein, pH, category), values in sorted(grouped.items()):
        mean, std, sem, n_values = group_stats(values)
        summaries.append({
            "protein": protein,
            "pH": pH,
            "category_key": category_key,
            "category": category,
            "mean_kcat_s": mean,
            "std_kcat_s": std,
            "sem_kcat_s": sem,
            "n_entries": n_values,
        })
    return summaries


def linear_series_summary(substrate_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    filtered = [
        row for row in substrate_rows
        if row["linear_series_name"] == "linear_primary_alcohols" and row["linear_series_index"] not in (None, "")
    ]
    grouped: dict[tuple[str, float], list[dict[str, object]]] = defaultdict(list)
    for row in filtered:
        grouped[(str(row["protein"]), float(row["pH"]))].append(row)
    summaries = []
    for (protein, pH), items in sorted(grouped.items()):
        items = sorted(items, key=lambda item: int(item["linear_series_index"]))
        corr, slope, intercept, n_points = safe_correlation(
            [item["carbon_number"] for item in items],
            [item["mean_kcat_s"] for item in items],
        )
        summaries.append({
            "protein": protein,
            "pH": pH,
            "series_name": "linear_primary_alcohols",
            "pearson_r_carbon_number": corr,
            "slope_kcat_per_carbon": slope,
            "intercept": intercept,
            "n_points": n_points,
        })
    return summaries


def write_substrate_templates(path: Path, descriptors: dict[str, SubstrateDescriptor], missing: set[str]) -> None:
    fieldnames = [
        "substrate", "substrate_class", "smiles", "carbon_number", "linear_series_name", "linear_series_index",
        "molecular_weight_da", "size_descriptor_name", "size_descriptor_value", "logp", "tpsa_a2",
        "rotatable_bonds", "hydroxyl_count", "alcohol_identity", "aromatic_ring_count", "acid_pka",
        "base_pka", "oxidation_potential_v", "oxidation_potential_source", "notes",
    ]
    rows = []
    for descriptor in sorted(descriptors.values(), key=lambda item: item.substrate):
        rows.append({
            "substrate": descriptor.substrate,
            "substrate_class": descriptor.substrate_class,
            "smiles": descriptor.smiles,
            "carbon_number": descriptor.carbon_number,
            "linear_series_name": descriptor.linear_series_name,
            "linear_series_index": finite_or_blank(descriptor.linear_series_index),
            "molecular_weight_da": descriptor.molecular_weight_da,
            "size_descriptor_name": descriptor.size_descriptor_name,
            "size_descriptor_value": descriptor.size_descriptor_value,
            "logp": descriptor.logp,
            "tpsa_a2": descriptor.tpsa_a2,
            "rotatable_bonds": descriptor.rotatable_bonds,
            "hydroxyl_count": descriptor.hydroxyl_count,
            "alcohol_identity": descriptor.alcohol_identity,
            "aromatic_ring_count": descriptor.aromatic_ring_count,
            "acid_pka": finite_or_blank(descriptor.acid_pka),
            "base_pka": finite_or_blank(descriptor.base_pka),
            "oxidation_potential_v": finite_or_blank(descriptor.oxidation_potential_v),
            "oxidation_potential_source": descriptor.oxidation_potential_source,
            "notes": descriptor.notes,
        })
    for substrate in sorted(missing):
        rows.append({
            "substrate": substrate,
            "substrate_class": "",
            "smiles": "",
            "carbon_number": "",
            "linear_series_name": "",
            "linear_series_index": "",
            "molecular_weight_da": "",
            "size_descriptor_name": "molecular_weight_da",
            "size_descriptor_value": "",
            "logp": "",
            "tpsa_a2": "",
            "rotatable_bonds": "",
            "hydroxyl_count": "",
            "alcohol_identity": "",
            "aromatic_ring_count": "",
            "acid_pka": "",
            "base_pka": "",
            "oxidation_potential_v": "",
            "oxidation_potential_source": "",
            "notes": "Fill this row for an unrecognized substrate",
        })
    write_dict_rows(path, rows, fieldnames)


def write_metal_templates(path: Path, descriptors: dict[str, MetalDescriptor], missing: set[str]) -> None:
    fieldnames = [
        "metal", "atomic_number", "charge", "ionic_radius_pm", "coordination_number",
        "f_electron_count", "is_lanthanide", "source", "notes",
    ]
    rows = []
    for descriptor in sorted(descriptors.values(), key=lambda item: item.metal):
        rows.append({
            "metal": descriptor.metal,
            "atomic_number": descriptor.atomic_number,
            "charge": descriptor.charge,
            "ionic_radius_pm": descriptor.ionic_radius_pm,
            "coordination_number": descriptor.coordination_number,
            "f_electron_count": descriptor.f_electron_count,
            "is_lanthanide": int(descriptor.is_lanthanide),
            "source": descriptor.source,
            "notes": descriptor.notes,
        })
    for metal in sorted(missing):
        rows.append({
            "metal": metal,
            "atomic_number": "",
            "charge": "",
            "ionic_radius_pm": "",
            "coordination_number": 6,
            "f_electron_count": "",
            "is_lanthanide": "",
            "source": "",
            "notes": "Fill this row for an unrecognized metal",
        })
    write_dict_rows(path, rows, fieldnames)


def write_metal_kinetics_template(path: Path, proteins: list[str], pH_values: list[float], metals: list[str]) -> None:
    fieldnames = ["protein", "pH", "metal", "apparent_halfmax_metal_uM", "apparent_metal_kd_uM", "source"]
    rows = []
    for protein in proteins:
        for pH in pH_values:
            for metal in metals:
                rows.append({
                    "protein": protein,
                    "pH": pH,
                    "metal": metal,
                    "apparent_halfmax_metal_uM": "",
                    "apparent_metal_kd_uM": "",
                    "source": "Requires metal titration data; not inferable from identity screen alone",
                })
    write_dict_rows(path, rows, fieldnames)


def write_analysis_notes(path: Path, proteins: list[str], pH_values: list[float], records: list[WellRecord]) -> None:
    lines = [
        "compare_proteins.py analysis notes",
        "",
        f"Proteins: {', '.join(proteins)}",
        f"pH values: {', '.join(f'{value:g}' for value in pH_values)}",
        f"Total non-control well rows: {len(records)}",
        "",
        "Built-in variables included in the output tables:",
        "- Substrates: class, carbon number, molecular weight, size descriptor, logP, logD at assay pH, predicted charge at assay pH, TPSA, rotatable bonds, hydroxyl count, alcohol identity, aromatic ring count, optional oxidation potential.",
        "- Metals: ionic radius (Shannon CN=6), atomic number, 4f electron count, charge, charge density, optional half-max metal concentration and apparent Kd.",
        "",
        "Variables not derivable from the provided screen alone:",
        "- Apparent metal Kd / half-max metal concentration require separate metal titration data.",
        "- Oxidation potential is left blank unless provided in a descriptor override CSV.",
        "",
        "Recommended next steps after the first run:",
        "1. Inspect *_condition_means.csv to confirm substrate and metal names matched the bundled descriptor table.",
        "2. Fill *_metal_kinetics_template.csv if you have titration-derived half-max or Kd values.",
        "3. Add oxidation potentials to *_substrate_descriptor_template.csv only if you have reference-quality values collected under a consistent convention.",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def plot_heatmaps(
    rows: list[ConditionMean],
    proteins: list[str],
    pH_values: list[float],
    substrate_descriptors: dict[str, SubstrateDescriptor],
    out_path: Path,
    vmax: float | None,
) -> None:
    plt.rcParams.update(NATURE_RC)
    substrates = sorted({row.substrate for row in rows}, key=lambda name: (resolve_substrate_descriptor(name, substrate_descriptors) is None, name))
    substrate_lookup = {name: resolve_substrate_descriptor(name, substrate_descriptors) for name in substrates}
    substrates = sorted(substrates, key=lambda name: (
        substrate_lookup[name].carbon_number if substrate_lookup[name] else 999,
        substrate_lookup[name].substrate_class if substrate_lookup[name] else "",
        name,
    ))
    metals = [metal for metal in METAL_ORDER if any(row.metal == metal for row in rows)]
    metals += sorted({row.metal for row in rows} - set(metals))
    lookup = {(row.protein, row.pH, row.substrate, row.metal): row.mean_kcat for row in rows}
    figure, axes = plt.subplots(len(pH_values), len(proteins), figsize=(9.2, max(3.0, 2.1 * len(pH_values))), squeeze=False)
    values = np.array([row.mean_kcat for row in rows], dtype=float)
    lower = min(0.0, float(np.nanmin(values)))
    upper = vmax if vmax is not None else float(np.nanmax(values))
    for row_index, pH in enumerate(pH_values):
        for col_index, protein in enumerate(proteins):
            axis = axes[row_index, col_index]
            grid = np.array([[lookup.get((protein, pH, substrate, metal), np.nan) for metal in metals] for substrate in substrates])
            image = axis.imshow(grid, aspect="auto", interpolation="nearest", cmap="viridis", vmin=lower, vmax=upper)
            axis.set_title(f"{protein}, pH {pH:g}")
            axis.set_xticks(range(len(metals)), labels=metals, rotation=90)
            axis.tick_params(axis="x", labelsize=7, length=0)
            axis.tick_params(axis="y", labelsize=7, length=0, pad=2)
            axis.set_yticks(range(len(substrates)), labels=substrates)
            for label in axis.get_yticklabels():
                label.set_horizontalalignment("right")
    figure.subplots_adjust(left=0.18, right=0.88, bottom=0.11, top=0.93, wspace=0.56, hspace=0.28)
    color_axis = figure.add_axes((0.90, 0.16, 0.02, 0.68))
    colorbar = figure.colorbar(image, cax=color_axis)
    colorbar.set_label("kcat (s$^{-1}$)")
    finalize_figure(figure, out_path)


def build_condition_legend_handles(
    proteins: list[str],
    pH_values: list[float],
    colors: list[tuple[float, float, float, float] | tuple[float, float, float]],
    marker_map: dict[str, str],
) -> list[Line2D]:
    handles: list[Line2D] = []
    for protein_index, protein in enumerate(proteins):
        color = colors[protein_index % len(colors)]
        for pH_index, pH in enumerate(pH_values):
            alpha = 0.35 + 0.65 * (pH_index / max(1, len(pH_values) - 1))
            handles.append(Line2D(
                [], [], linestyle="none", marker=marker_map[protein], markersize=5,
                markerfacecolor=color, markeredgewidth=0, alpha=alpha,
                label=f"{protein}, pH {pH:g}",
            ))
    return handles


def add_stacked_condition_legends(
    figure,
    proteins: list[str],
    pH_values: list[float],
    colors: list[tuple[float, float, float, float] | tuple[float, float, float]],
    marker_map: dict[str, str],
) -> None:
    for legend_index, protein in enumerate(proteins):
        handles: list[Line2D] = []
        protein_index = proteins.index(protein)
        color = colors[protein_index % len(colors)]
        for pH_index, pH in enumerate(pH_values):
            alpha = 0.35 + 0.65 * (pH_index / max(1, len(pH_values) - 1))
            handles.append(Line2D(
                [], [], linestyle="none", marker=marker_map[protein], markersize=5,
                markerfacecolor=color, markeredgewidth=0, alpha=alpha,
                label=f"{protein} pH {pH:g}",
            ))
        legend = figure.legend(
            handles=handles,
            frameon=False,
            loc="upper center",
            ncol=max(1, len(pH_values)),
            bbox_to_anchor=(0.5, 0.995 - legend_index * 0.045),
            handletextpad=0.4,
            columnspacing=1.0,
        )
        figure.add_artist(legend)


def plot_descriptor_panels(
    rows: list[dict[str, object]],
    proteins: list[str],
    pH_values: list[float],
    descriptors: list[str],
    xlabels: dict[str, str],
    out_path: Path,
) -> None:
    plt.rcParams.update(NATURE_RC)
    colors = plt.get_cmap("tab10").colors
    figure, axes = plt.subplots(int(math.ceil(len(descriptors) / 2)), 2, figsize=(7.1, 2.4 * int(math.ceil(len(descriptors) / 2))), squeeze=False)
    marker_map = {protein: marker for protein, marker in zip(proteins, ["o", "s", "^", "D", "v", "P"])}
    flat_axes = axes.ravel()
    for axis, descriptor in zip(flat_axes, descriptors):
        pretty_label = xlabels.get(descriptor, descriptor.replace("_", " "))
        for protein_index, protein in enumerate(proteins):
            for pH_index, pH in enumerate(pH_values):
                subset = [row for row in rows if row["protein"] == protein and float(row["pH"]) == pH]
                x_array, y_array = finite_pairs([item[descriptor] for item in subset], [item["mean_kcat_s"] for item in subset])
                if x_array.size == 0:
                    continue
                color = colors[protein_index % len(colors)]
                alpha = 0.35 + 0.65 * (pH_index / max(1, len(pH_values) - 1))
                axis.scatter(x_array, y_array, s=16, marker=marker_map[protein], color=color, alpha=alpha,
                             label=f"{protein}, pH {pH:g}")
        axis.set_xlabel("")
        axis.set_ylabel("mean kcat (s$^{-1}$)")
        axis.set_title(pretty_label)
    for axis in flat_axes[len(descriptors):]:
        axis.axis("off")
    add_stacked_condition_legends(figure, proteins, pH_values, list(colors), marker_map)
    figure.tight_layout(rect=(0, 0, 1, 0.89))
    finalize_figure(figure, out_path)


def plot_protein_scatter(pair_rows: list[dict[str, object]], proteins: list[str], out_path: Path) -> None:
    if len(proteins) != 2 or not pair_rows:
        return
    left, right = proteins
    x_key = f"{left}_kcat_s"
    y_key = f"{right}_kcat_s"
    plt.rcParams.update(NATURE_RC)
    figure, axis = plt.subplots(figsize=(4.2, 4.0))
    cmap = plt.get_cmap("viridis")
    pH_values = sorted({float(row["pH"]) for row in pair_rows})
    pH_min = min(pH_values)
    pH_span = max(pH_values) - pH_min if len(pH_values) > 1 else 1.0
    legend_handles: list[Line2D] = []
    for pH in pH_values:
        color = cmap((pH - pH_min) / pH_span)
        legend_handles.append(Line2D([], [], linestyle="none", marker="o", markersize=5,
                                     markerfacecolor=color, markeredgewidth=0, alpha=0.8,
                                     label=f"pH {pH:g}"))
    for row in pair_rows:
        pH = float(row["pH"])
        color = cmap((pH - pH_min) / pH_span)
        axis.scatter(float(row[x_key]), float(row[y_key]), s=14, color=color, alpha=0.7)
    finite_x, finite_y = finite_pairs([row[x_key] for row in pair_rows], [row[y_key] for row in pair_rows])
    if finite_x.size:
        low = min(float(np.min(finite_x)), float(np.min(finite_y)))
        high = max(float(np.max(finite_x)), float(np.max(finite_y)))
        padding = 0.04 * (high - low if high > low else max(high, 1.0))
        axis.set_xlim(low - padding, high + padding)
        axis.set_ylim(low - padding, high + padding)
        axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(f"{left} mean kcat (s$^{-1}$)")
    axis.set_ylabel(f"{right} mean kcat (s$^{-1}$)")
    axis.set_title(f"{right} vs {left} matched-condition comparison")
    axis.legend(handles=legend_handles, frameon=False, loc="best")
    figure.tight_layout()
    finalize_figure(figure, out_path)


def finalize_figure(figure, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = out_path.resolve()
    png_path = out_path.with_suffix(".png").resolve()
    figure.savefig(str(pdf_path), format="pdf")
    try:
        with tempfile.NamedTemporaryFile(suffix=".png", dir=png_path.parent, delete=False) as handle:
            temp_png = Path(handle.name)
        try:
            figure.savefig(str(temp_png), format="png", dpi=600)
            if os.name == "nt":
                os.replace(str(temp_png), "\\\\?\\" + str(png_path))
            else:
                os.replace(str(temp_png), str(png_path))
        finally:
            if temp_png.exists():
                temp_png.unlink()
    except OSError as exc:
        print(f"WARNING: could not write PNG {png_path.name}: {exc}")
    plt.close(figure)


def main() -> None:
    parser = build_argparser()
    args = parser.parse_args()
    if args.kcat_upper is not None and (not np.isfinite(args.kcat_upper) or args.kcat_upper <= 0):
        parser.error("--kcat-upper must be finite and greater than 0.")

    conditions = resolve_conditions(args, parser)
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    records: list[WellRecord] = []
    for protein, pH, csv_path in conditions:
        records.extend(read_wells_csv(csv_path, protein, pH))
    if not records:
        raise SystemExit("No usable non-control kcat rows were found in the supplied files.")

    proteins = sorted({record.protein for record in records})
    pH_values = sorted({record.pH for record in records})
    metals_present = sorted({record.metal for record in records}, key=lambda metal: (METAL_ORDER.index(metal) if metal in METAL_ORDER else 999, metal))
    substrates_present = sorted({record.substrate for record in records})

    substrate_descriptors = load_substrate_descriptors(args.substrate_descriptors)
    metal_descriptors = load_metal_descriptors(args.metal_descriptors)
    metal_kinetics = load_metal_kinetics(args.metal_kinetics)

    missing_substrates = {substrate for substrate in substrates_present if resolve_substrate_descriptor(substrate, substrate_descriptors) is None}
    missing_metals = {metal for metal in metals_present if metal not in metal_descriptors}
    write_substrate_templates(outdir / f"{args.label}_substrate_descriptor_template.csv", substrate_descriptors, missing_substrates)
    write_metal_templates(outdir / f"{args.label}_metal_descriptor_template.csv", metal_descriptors, missing_metals)
    write_metal_kinetics_template(outdir / f"{args.label}_metal_kinetics_template.csv", proteins, pH_values, metals_present)
    if missing_substrates or missing_metals:
        problems = []
        if missing_substrates:
            problems.append(f"Unrecognized substrates: {sorted(missing_substrates)}")
        if missing_metals:
            problems.append(f"Unrecognized metals: {sorted(missing_metals)}")
        problems.append("Fill the exported template CSVs and rerun with --substrate-descriptors and/or --metal-descriptors.")
        raise SystemExit(" ".join(problems))

    condition_rows = aggregate_condition_means(records, substrate_descriptors, metal_descriptors, metal_kinetics)
    substrate_rows = aggregate_by_substrate(condition_rows)
    metal_rows = aggregate_by_metal(condition_rows)
    pair_rows = matched_protein_rows(condition_rows)

    write_wells_long(outdir / f"{args.label}_wells_long.csv", records, substrate_descriptors, metal_descriptors)
    write_condition_means(outdir / f"{args.label}_condition_means.csv", condition_rows)
    write_dict_rows(
        outdir / f"{args.label}_substrate_summary.csv",
        substrate_rows,
        [
            "protein", "pH", "substrate", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_metals",
            "substrate_class", "carbon_number", "molecular_weight_da", "size_descriptor_name",
            "size_descriptor_value", "logp", "predicted_charge", "logd", "tpsa_a2", "rotatable_bonds",
            "hydroxyl_count", "alcohol_identity", "aromatic_ring_count", "oxidation_potential_v",
            "linear_series_name", "linear_series_index",
        ],
    )
    write_dict_rows(
        outdir / f"{args.label}_metal_summary.csv",
        metal_rows,
        [
            "protein", "pH", "metal", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_substrates",
            "metal_charge", "ionic_radius_pm", "coordination_number", "atomic_number", "f_electron_count",
            "charge_density_charge_per_a2", "is_lanthanide", "apparent_halfmax_metal_uM",
            "apparent_metal_kd_uM", "metal_kinetics_source",
        ],
    )
    if pair_rows:
        write_dict_rows(outdir / f"{args.label}_protein_pairs.csv", pair_rows, list(pair_rows[0].keys()))

    substrate_descriptor_names = [
        "size_descriptor_value", "molecular_weight_da", "carbon_number", "logp", "logd", "tpsa_a2",
        "rotatable_bonds", "hydroxyl_count", "aromatic_ring_count", "predicted_charge", "oxidation_potential_v",
    ]
    metal_descriptor_names = [
        "ionic_radius_pm", "atomic_number", "f_electron_count", "charge_density_charge_per_a2",
        "apparent_halfmax_metal_uM", "apparent_metal_kd_uM",
    ]
    substrate_correlations = descriptor_correlations(substrate_rows, "substrate", substrate_descriptor_names)
    substrate_correlations += subset_correlations(
        substrate_rows,
        "linear_primary_alcohols",
        lambda row: row["linear_series_name"] == "linear_primary_alcohols",
        ["carbon_number", "molecular_weight_da", "logp", "logd"],
    )
    substrate_correlations += subset_correlations(
        substrate_rows,
        "alcohol_and_aldehyde_only",
        lambda row: row["substrate_class"] in {"alcohol", "aldehyde", "amino_alcohol"},
        ["molecular_weight_da", "carbon_number", "logp", "logd", "tpsa_a2", "rotatable_bonds", "predicted_charge"],
    )
    metal_correlations = descriptor_correlations(metal_rows, "metal", metal_descriptor_names)
    metal_correlations += subset_correlations(
        metal_rows,
        "lanthanides_only",
        lambda row: bool(row["is_lanthanide"]),
        ["ionic_radius_pm", "atomic_number", "f_electron_count", "charge_density_charge_per_a2", "apparent_halfmax_metal_uM", "apparent_metal_kd_uM"],
    )
    write_dict_rows(
        outdir / f"{args.label}_substrate_correlations.csv",
        substrate_correlations,
        ["protein", "pH", "axis", "descriptor", "subset", "pearson_r", "slope", "intercept", "n_points"],
    )
    write_dict_rows(
        outdir / f"{args.label}_metal_correlations.csv",
        metal_correlations,
        ["protein", "pH", "axis", "descriptor", "subset", "pearson_r", "slope", "intercept", "n_points"],
    )
    write_dict_rows(
        outdir / f"{args.label}_substrate_class_summary.csv",
        categorical_group_summary(substrate_rows, "substrate_class"),
        ["protein", "pH", "category_key", "category", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_entries"],
    )
    write_dict_rows(
        outdir / f"{args.label}_alcohol_identity_summary.csv",
        categorical_group_summary(substrate_rows, "alcohol_identity"),
        ["protein", "pH", "category_key", "category", "mean_kcat_s", "std_kcat_s", "sem_kcat_s", "n_entries"],
    )
    write_dict_rows(
        outdir / f"{args.label}_linear_series_summary.csv",
        linear_series_summary(substrate_rows),
        ["protein", "pH", "series_name", "pearson_r_carbon_number", "slope_kcat_per_carbon", "intercept", "n_points"],
    )
    write_analysis_notes(outdir / f"{args.label}_analysis_notes.txt", proteins, pH_values, records)

    if not args.no_plots:
        plot_heatmaps(
            condition_rows,
            proteins,
            pH_values,
            substrate_descriptors,
            outdir / f"{args.label}_heatmaps.pdf",
            args.kcat_upper,
        )
        plot_descriptor_panels(
            substrate_rows,
            proteins,
            pH_values,
            ["molecular_weight_da", "carbon_number", "logp", "logd", "tpsa_a2", "rotatable_bonds", "hydroxyl_count", "predicted_charge"],
            {
                "molecular_weight_da": "substrate MW (Da)",
                "carbon_number": "substrate carbon number",
                "logp": "substrate logP",
                "logd": "substrate logD at assay pH",
                "tpsa_a2": "substrate TPSA (Å$^2$)",
                "rotatable_bonds": "rotatable bond count",
                "hydroxyl_count": "hydroxyl count",
                "predicted_charge": "predicted charge at assay pH",
            },
            outdir / f"{args.label}_substrate_descriptors.pdf",
        )
        plot_descriptor_panels(
            metal_rows,
            proteins,
            pH_values,
            ["ionic_radius_pm", "atomic_number", "f_electron_count", "charge_density_charge_per_a2"],
            {
                "ionic_radius_pm": "ionic radius (pm)",
                "atomic_number": "atomic number",
                "f_electron_count": "4f electron count",
                "charge_density_charge_per_a2": "charge density (charge Å$^{-2}$)",
            },
            outdir / f"{args.label}_metal_descriptors.pdf",
        )
        if pair_rows:
            plot_protein_scatter(pair_rows, proteins, outdir / f"{args.label}_protein_scatter.pdf")

    print(f"Wrote outputs to {outdir}")
    print(f"Proteins: {proteins}")
    print(f"pH values: {[f'{value:g}' for value in pH_values]}")
    print(f"Metals: {len(metals_present)}")
    print(f"Substrates: {len(substrates_present)}")


if __name__ == "__main__":
    main()