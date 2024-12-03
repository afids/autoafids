from __future__ import annotations

import csv
from os import PathLike
from pathlib import Path

from numpy.typing import NDArray

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
    Path(__file__).parent / ".." / ".." / "resources" / "tpl-MNI152NLin2009cAsym_res-01_T1w.fcsv"
)


def afids_to_fcsv(
    afid_coords: dict[int, NDArray],
    fcsv_output: PathLike[str] | str,
) -> None:
    """
    AFIDS to Slicer-compatible .fcsv file.
    
    Parameters
    ----------

        afids_coords :: dict
            AFIDS coordinates

        fcsv_output :: str
            Path to output fcsv file

    Returns
    -------
        None
    """
    # Read in fcsv template
    with FCSV_TEMPLATE.open(encoding="utf-8", newline="") as fcsv_file:
        header = [fcsv_file.readline() for _ in range(3)]
        reader = csv.DictReader(fcsv_file, fieldnames=AFIDS_FIELDNAMES)
        fcsv = list(reader)

    # Loop over fiducials
    for idx, row in enumerate(fcsv):
        # Update fcsv, skipping header
        label = idx + 1
        row["x"] = afid_coords[label][0]
        row["y"] = afid_coords[label][1]
        row["z"] = afid_coords[label][2]

    # Write output fcsv
    with Path(fcsv_output).open("w", encoding="utf-8", newline="") as out_fcsv_file:
        for line in header:
            out_fcsv_file.write(line)
        writer = csv.DictWriter(out_fcsv_file, fieldnames=AFIDS_FIELDNAMES)
        for row in fcsv:
            writer.writerow(row)