"""nnlm_to_fcsv.py
Snakemake script: convert nnLM voxel-space JSON output to a Slicer
FCSV file with RAS world coordinates (mm).

nnLM JSON convention
--------------------
  {"1": {"coordinates": [x, y, z], "likelihood": 0.87}, ...}
  where x = column (SimpleITK axis 0), y = row, z = slice.

NiBabel affine convention
--------------------------
  affine @ [i, j, k, 1]  with  i = slice-index, j = row, k = col
  → maps (z_sitk, y_sitk, x_sitk) to RAS world coords.
"""

import csv
import json
from pathlib import Path

import nibabel as nib
import numpy as np

# ── Load nnLM voxel-space coordinates ────────────────────────────────────────
with open(snakemake.input.coords_json) as fh:
    nnlm_raw = json.load(fh)   # {"1": {"coordinates": [x,y,z], ...}, ...}

# ── Load image affine (nibabel i,j,k = slice,row,col) ────────────────────────
img = nib.load(snakemake.input.t1w)
affine = img.affine   # 4×4 voxel→RAS matrix

# ── Convert voxel → RAS world coords ─────────────────────────────────────────
# SimpleITK (x,y,z) == (col, row, slice)  →  nibabel (i,j,k) == (slice, row, col)
afid_world: dict[int, np.ndarray] = {}
for afid_str, entry in nnlm_raw.items():
    sx, sy, sz = entry["coordinates"]            # SimpleITK x,y,z
    # NiBabel expects (x, y, z) index mapping identical to SimpleITK ordering
    nib_ijk = np.array([sx, sy, sz, 1.0])
    # 1) Get RAS world coordinates from the NIfTI affine (this matches true FCSV space)
    world = (affine @ nib_ijk)[:3]
    afid_world[int(afid_str)] = world

# ── Read FCSV template (header + 32 landmark rows) ───────────────────────────
fcsv_template = Path(snakemake.params.fcsv_template)

FIELDNAMES = [
    "id", "x", "y", "z",
    "ow", "ox", "oy", "oz",
    "vis", "sel", "lock",
    "label", "desc", "associatedNodeID",
]

with fcsv_template.open(encoding="utf-8", newline="") as fh:
    # Preserve the 3-line Slicer header verbatim
    header_lines = [fh.readline() for _ in range(3)]
    reader = csv.DictReader(fh, fieldnames=FIELDNAMES)
    rows = list(reader)

# ── Fill in predicted coords ──────────────────────────────────────────────────
for idx, row in enumerate(rows):
    afid_label = idx + 1
    coords = afid_world.get(afid_label, np.zeros(3))
    row["x"] = f"{coords[0]:.6f}"
    row["y"] = f"{coords[1]:.6f}"
    row["z"] = f"{coords[2]:.6f}"

# ── Write output FCSV ─────────────────────────────────────────────────────────
out_path = Path(snakemake.output.fcsv)
out_path.parent.mkdir(parents=True, exist_ok=True)

with out_path.open("w", encoding="utf-8", newline="") as fh:
    for line in header_lines:
        fh.write(line)
    writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, extrasaction="ignore")
    for row in rows:
        writer.writerow(row)
