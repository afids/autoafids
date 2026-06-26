# ruff: noqa
"""
gather_afids.py
===============

Snakemake script: collect all 32 per-AFID coord files and write the
combined FCSV.

Called by the ``applyfidmodel_gather`` rule after all 32
``applyfidmodel_single`` jobs have completed.

Snakemake I/O
-------------
    input:
        coords – list of 32 plain-text coord files (x y z), one per AFID
    output:
        fcsv   – combined Slicer-compatible FCSV with all 32 AFIDs
"""

import csv
import warnings
from pathlib import Path
from typing import Dict

import numpy as np

warnings.filterwarnings("ignore")


# ===================================================================
# FCSV writer  (identical to apply_with_prior.py)
# ===================================================================

AFIDS_FIELDNAMES = [
    "id",
    "x",
    "y",
    "z",
    "ow",
    "ox",
    "oy",
    "oz",
    "vis",
    "sel",
    "lock",
    "label",
    "desc",
    "associatedNodeID",
]

FCSV_TEMPLATE = (
    Path(__file__).parent
    / ".."
    / ".."
    / "resources"
    / "tpl-MNI152NLin2009cAsym_res-01_T1w.fcsv"
)


def afids_to_fcsv(afid_coords: Dict[int, np.ndarray], fcsv_output) -> None:
    with FCSV_TEMPLATE.open(encoding="utf-8", newline="") as fcsv_file:
        header = [fcsv_file.readline() for _ in range(3)]
        reader = csv.DictReader(fcsv_file, fieldnames=AFIDS_FIELDNAMES)
        fcsv = list(reader)

    for idx, row in enumerate(fcsv):
        label = idx + 1
        c = afid_coords.get(label, np.array([0.0, 0.0, 0.0]))
        row["x"] = c[0]
        row["y"] = c[1]
        row["z"] = c[2]

    with Path(fcsv_output).open("w", encoding="utf-8", newline="") as out:
        for line in header:
            out.write(line)
        writer = csv.DictWriter(out, fieldnames=AFIDS_FIELDNAMES)
        for row in fcsv:
            writer.writerow(row)


# ===================================================================
# Snakemake entry point
# ===================================================================

# input.coords is a list of coord files sorted by afid wildcard order.
# We infer the AFID number from the filename (afid-XX in the path).
afid_coords: Dict[int, np.ndarray] = {}

for coord_path in snakemake.input.coords:  # noqa: F821
    path = Path(coord_path)
    # Extract afid number from filename, e.g. "sub-001_afid-03_coord.txt" -> 3
    for part in path.stem.split("_"):
        if part.startswith("afid-"):
            fid = int(part.split("-")[1])
            break
    else:
        raise ValueError(f"Cannot parse AFID number from filename: {path.name}")

    with open(coord_path) as f:
        x, y, z = map(float, f.read().strip().split())
    afid_coords[fid] = np.array([x, y, z])
    print(f"  [gather] AFID {fid:02d}  world=[{x:.1f}, {y:.1f}, {z:.1f}]")

print(f"\n  Writing combined FCSV with {len(afid_coords)} AFIDs...")
Path(snakemake.output.fcsv).parent.mkdir(parents=True, exist_ok=True)  # noqa: F821
afids_to_fcsv(afid_coords, snakemake.output.fcsv)  # noqa: F821
print(f"  Done → {snakemake.output.fcsv}")  # noqa: F821
