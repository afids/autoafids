from __future__ import annotations

import base64
import csv
import os
import random as rand
import re
from io import BytesIO
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from appdirs import AppDirs
from joblib import Parallel, delayed
from numpy.typing import NDArray

# Dictionary for AFID labels
afids_labels = {
    1: "AC",
    2: "PC",
    3: "ICS",
    4: "PMJ",
    5: "SIPF",
    6: "RSLMS",
    7: "LSLMS",
    8: "RILMS",
    9: "LILMS",
    10: "CUL",
    11: "IMS",
    12: "RMB",
    13: "LMB",
    14: "PG",
    15: "RLVAC",
    16: "LLVAC",
    17: "RLVPC",
    18: "LLVPC",
    19: "GENU",
    20: "SPLE",
    21: "RALTH",
    22: "LALTH",
    23: "RSAMTH",
    24: "LSAMTH",
    25: "RIAMTH",
    26: "LIAMTH",
    27: "RIGO",
    28: "LIGO",
    29: "RVOH",
    30: "LVOH",
    31: "ROSF",
    32: "LOSF",
}


def fcsvtodf(fcsv_path):
    """
    Convert a .fcsv file (assumes RAS coordinate system)
    to a ML-friendly dataframe and return the cleaned xyz coordinates.

    Parameters:
    - fcsv_path: str, path to the .fcsv file

    Returns:
    - df_xyz_clean: pandas.DataFrame with cleaned x, y, z coordinates
    - num_points: int, number of fiducial points
    """

    # Extract the subject ID from the file path (naming is in bids-like)
    subject_id = re.search(r"(sub-\w+)", fcsv_path).group(1)

    # Read in .fcsv file, skip header
    df_raw = pd.read_table(fcsv_path, sep=",", header=2)

    # Extract the x, y, z coordiantes and store them in data science
    # friendly format (i.e., features in cols and subject in rows)
    df_xyz = df_raw[["x", "y", "z"]].melt().transpose()

    # Use number of row in fcsv to make number points
    colnames = [
        f"{axis}_{i % int(df_raw.shape[0]) + 1}"
        for axis in ["x", "y", "z"]
        for i in range(int(df_raw.shape[0]))
    ]

    # Reassign features to be descriptive of coordinate
    df_xyz.columns = colnames

    # clean dataframe and pin correct subject name
    df_xyz_clean = df_xyz.drop("variable", axis=0)
    df_xyz_clean = df_xyz_clean.rename(index={"value": subject_id})
    df_xyz_clean = df_xyz_clean.astype(float)

    return df_xyz_clean, df_raw.shape[0]


def dftodfml(fcsvdf):
    """
    Convert a datafrane (assumes RAS coordinate system)
    to a ML-friendly dataframe and return the cleaned xyz coordinates.

    Parameters:
    - fcsvdf: pandas.DataFrame

    Returns:
    - df_xyz_clean: pandas.DataFrame with cleaned x, y, z coordinates
    - num_points: int, number of fiducial points
    """

    # Extract the x, y, z coordiantes and store them
    # in data science friendly format
    # (i.e., features in cols and subject in rows)
    df_xyz = fcsvdf[["x", "y", "z"]].melt().transpose()

    # Use number of row in fcsv to make number points
    colnames = [
        f"{axis}_{i % int(fcsvdf.shape[0]) + 1}"
        for axis in ["x", "y", "z"]
        for i in range(int(fcsvdf.shape[0]))
    ]

    # Reassign features to be descriptive of coordinate
    df_xyz.columns = colnames

    # clean dataframe and pin correct subject name
    df_xyz_clean = df_xyz.drop("variable", axis=0)
    df_xyz_clean = df_xyz_clean.astype(float)

    return df_xyz_clean, fcsvdf.shape[0]


def get_fiducial_index(fid):
    """
    Retrieve the index corresponding to the fiducial name or integer.

    Parameters:
    - fid: str or int, fiducial identifier (name or index)

    Returns:
    - int, corresponding fiducial index
    """

    if isinstance(fid, str):
        for idx, name in afids_labels.items():
            if name == fid:
                return idx
    elif isinstance(fid, int):
        return fid
    raise ValueError("Invalid fiducial identifier.")


def compute_distance(fcsv_path, fid1, fid2):
    """
    Compute the Euclidean distance between two fiducials.

    Parameters:
    - fcsv_path: str, path to the .fcsv file
    - fid1, fid2: str or int, fiducial identifiers

    Returns:
    - xyz_diff: numpy.array, difference in x, y, z coordinates
    - distance: float, Euclidean distance between fiducials
    """

    # Retrieve indices of the fiducials
    index1, index2 = get_fiducial_index(fid1), get_fiducial_index(fid2)

    # Load dataframe from the fcsv file
    df = fcsvtodf(fcsv_path)[0]

    # Extract x, y, z coordinates into numpy arrays
    coords1 = df[[f"x_{index1}", f"y_{index1}", f"z_{index1}"]].to_numpy()
    coords2 = df[[f"x_{index2}", f"y_{index2}", f"z_{index2}"]].to_numpy()

    # Compute the difference as a numpy array
    xyz_diff = coords1 - coords2

    # Compute the Euclidean distance
    distance = np.linalg.norm(xyz_diff)

    return xyz_diff.flatten(), distance


def compute_average(fcsv_path, fid1, fid2):
    """
    Compute the average position between two fiducials.

    Parameters:
    - fcsv_path: str, path to the .fcsv file
    - fid1, fid2: str or int, fiducial identifiers

    Returns:
    - xyz_average: numpy.array, average coordinates (x, y, z) between fiducials
    """

    # Retrieve indices of the fiducials
    index1, index2 = get_fiducial_index(fid1), get_fiducial_index(fid2)

    # Load dataframe from the fcsv file
    df = fcsvtodf(fcsv_path)[0]

    # Extract x, y, z coordinates into numpy arrays
    coords1 = df[[f"x_{index1}", f"y_{index1}", f"z_{index1}"]].to_numpy()
    coords2 = df[[f"x_{index2}", f"y_{index2}", f"z_{index2}"]].to_numpy()

    # Compute the average as a numpy array
    xyz_average = (coords1 + coords2) / 2

    return xyz_average.flatten()


def generate_slicer_file(matrix, filename):
    """
    Generate a .txt transformation file for 3D Slicer from a 4x4 matrix.

    Parameters:
    - matrix: np.ndarray, 4x4 transformation matrix
    - output_path: str, path to store .txt file
    """
    d = np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    ras_inmatrix = np.linalg.inv(matrix)

    lps_inmatrix = d @ ras_inmatrix @ d

    # Extract rotation/scale and translation components
    rotation_scale = lps_inmatrix[0:3, 0:3].flatten()
    translation = lps_inmatrix[0:3, 3]

    # Format the content of the .tfm file
    tfm_content = "#Insight Transform File V1.0\n"
    tfm_content += "#Transform 0\n"
    tfm_content += "Transform: AffineTransform_double_3_3\n"
    tfm_content += (
        "Parameters: "
        + " ".join(map(str, rotation_scale))
        + " "
        + " ".join(map(str, translation))
        + "\n"
    )
    tfm_content += "FixedParameters: 0 0 0\n"

    # Write the content to the specified file
    with open(filename, "w") as file:
        file.write(tfm_content)


def acpcmatrix(
    fcsv_path,
    midline,
    center_on_mcp=False,
    write_matrix=True,
    transform_file_name=None
):
    """
    Computes a 4x4 transformation matrix aligning with the AC-PC axis.

    Parameters:
    - fcsv_path: str
        path to the .fcsv file
    - midline: str or int
        one of the 10 midline AFID points.
    - center_on_mcp: bool
        If True, adds translation element
        to the ACPC matrix (centering on MCP).

    Returns:
    - matrix: np.ndarray, A 4x4 affine transformation matrix.
    """

    # A-P axis
    acpc = compute_distance(fcsv_path, "AC", "PC")  # vector from PC to AC
    # unit vector defining anterior and posterior axis
    y_axis = acpc[0] / acpc[1]

    # R-L axis
    lataxis = compute_distance(fcsv_path, midline, "AC")[
        0
    ]  # vector from AC to midline point
    x_axis = np.cross(y_axis, lataxis)  # vector defining left and right axis
    # unit vector defining left and right axis
    x_axis /= np.linalg.norm(x_axis)

    # S-I axis
    z_axis = np.cross(x_axis, y_axis)
    z_axis /= np.linalg.norm(z_axis)

    # Rotation matrix
    rotation = np.vstack([x_axis, y_axis, z_axis])
    translation = np.array([0, 0, 0])

    # Build 4x4 matrix
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation

    # Orientation correction for matrix
    # is midpoint placed below ACPC (i.e., at PMJ)
    if matrix[0][0] < 0:
        matrix = matrix * np.array(
            [[-1, -1, -1, 0], [1, 1, 1, 0], [-1, -1, -1, 0], [0, 0, 0, 1]]
        )

    if center_on_mcp:  # need to compute MCP AFTER rotation is applied
        mcp = compute_average(fcsv_path, "AC", "PC")  # MCP in native
        matrix[:3, 3] = -np.dot(
            matrix[:3, :3], mcp
        )  # MCP after rotation; negative because we set MCP to (0,0,0)

    if write_matrix:
        generate_slicer_file(matrix, transform_file_name)

    return matrix


def mcp_origin(df_afids):
    """
    sets MCP as the origin point for each of the subjects

    Parameters:
    - df_afids: pandas.DataFrame

    Returns:
    - df_afids_ori_mcp: pandas.DataFrame, mcp centered points.
    - mcp coordinates: tuple, x,y,z coordiantes of mcp point.
    """
    # extract MCP coordinates; defined as average point between AC and PC
    mcp_x = (df_afids["x_1"] + df_afids["x_2"]) / 2
    mcp_y = (df_afids["y_1"] + df_afids["y_2"]) / 2
    mcp_z = (df_afids["z_1"] + df_afids["z_2"]) / 2

    # subtract MCP coordinates from afids at appropriate coords
    df_afids_mcpx = df_afids.transpose()[0:32] - mcp_x
    df_afids_mcpy = df_afids.transpose()[32:64] - mcp_y
    df_afids_mcpz = df_afids.transpose()[64:98] - mcp_z

    # concat the three coords and take transpose
    frames = [df_afids_mcpx, df_afids_mcpy, df_afids_mcpz]
    df_afids_mcp = pd.concat(frames)
    df_afids_ori_mcp = df_afids_mcp.transpose()
    df_afids_ori_mcp = df_afids_ori_mcp.astype(float)

    return df_afids_ori_mcp, (
        np.array([mcp_x.to_numpy(), mcp_y.to_numpy(), mcp_z.to_numpy()])
    )


def make_zero(num, threshold=0.0001):
    if abs(num) < threshold:
        return 0
    else:
        return num


def transform_afids(fcsv_path, slicer_tfm, midpoint):
    """
    Computes and applies an AC-PC transformation to
    AFID files and saves an fcsv with transfromed coordinates.

    Parameters:
    - fcsv_path: str
        path to the .fcsv file for coordinates to be transformed
    - midline_fcsv: str
        path to the .fcsv file which has the midline coordinates
    - output_dir: str
        path of the directory to store transformed
        .fcsv file (if not specified, no output fcsv will be written)
    - midpoint: str or int, any midline AFID point
    - slicer_tfm: str, path to .txt file for slicer ACPC transform

    Returns:
    - tcoords: np.ndarray, transformed coordinates
    """

    # Compute the 4x4 AC-PC transformation matrix
    xfm_txt = acpcmatrix(
        fcsv_path=fcsv_path, midline=midpoint, transform_file_name=slicer_tfm
    )

    # Read coordinates from the file
    fcsv_df = pd.read_table(fcsv_path, sep=",", header=2)

    # Copy coordinates and apply transformation
    pretfm = fcsv_df.loc[:, ["x", "y", "z"]].copy()
    pretfm["coord"] = 1  # Add a fourth column for homogeneous coordinates
    coords = pretfm.to_numpy()
    thcoords = xfm_txt @ coords.T
    tcoords = thcoords.T[:, :3]  # Remove homogeneous coordinates

    # Substitute new coordinates
    if tcoords.shape == (len(fcsv_df), 3):
        fcsv_df.loc[:, "x"] = tcoords[:, 0]
        fcsv_df.loc[:, "y"] = tcoords[:, 1]
        fcsv_df.loc[:, "z"] = tcoords[:, 2]
    else:
        raise ValueError(
            "New coordinates do not match " +
            "the number of rows in the original fcsv."
        )

    return fcsv_df, xfm_txt


def save_single_slice_in_memory(
        data, x, y, z,
        offset,
        zoom_radius,
        show_crosshairs
    ):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for i, (slice_data, coord, title) in enumerate(
        zip(
            [data[x + offset, :, :],
            data[:, y + offset, :],
            data[:, :, z + offset]],
            [(y, z), (x, z), (x, y)],
            ["Sagittal", "Coronal", "Axial"],
        )
    ):
        axes[i].imshow(slice_data.T, origin="lower", cmap="gray")
        if offset == 0 and show_crosshairs:
            axes[i].axhline(y=coord[1], color="r", lw=1)
            axes[i].axvline(x=coord[0], color="r", lw=1)
        axes[i].set_xlim(coord[0] - zoom_radius, coord[0] + zoom_radius)
        axes[i].set_ylim(coord[1] - zoom_radius, coord[1] + zoom_radius)
        axes[i].set_title(title)
        axes[i].axis("off")

    buffer = BytesIO()
    plt.savefig(buffer, format="png", bbox_inches="tight")
    plt.close(fig)
    buffer.seek(0)

    return base64.b64encode(buffer.read()).decode()


def save_mri_slices_as_images(
        data, x, y, z,
        jitter,
        zoom_radius,
        show_crosshairs=True):
    return Parallel(n_jobs=-1)(
        delayed(save_single_slice_in_memory)(
            data, x, y, z, offset, zoom_radius, show_crosshairs
        )
        for offset in range(-jitter, jitter + 1)
    )


def extract_coordinates_from_fcsv(file_path, label_description):
    df = pd.read_csv(file_path, comment="#", header=None)
    row = df[df[11] == label_description]
    return tuple(row.iloc[0, 1:4]) if not row.empty else None


def generate_html_with_keypress(
    subject_images, reference_images, output_html="mri_viewer.html"
):
    """Generate an interactive HTML viewer with sticky instructions."""
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>MRI Viewer</title>
        <style>
            body { display: flex;
                   font-family: Arial,
                   sans-serif;
                   margin: 0;
                   padding: 0; }
            .instructions {
                width: 20%; padding: 20px; background-color: #f4f4f4;
                border-right: 2px solid #ddd;
                position: sticky; top: 0; height: 100vh; overflow-y: auto;
                box-sizing: border-box;
            }
            .instructions h2 { margin-top: 0; }
            .viewer { flex: 1; padding: 20px; text-align: center; }
            .slider { width: 80%; margin: 20px auto; }
            .image { display:
                     block;
                     margin: 0 auto;
                     max-width: 80%;
                     transition: opacity 0.2s ease-in-out; }
        </style>
    </head>
    <body>
        <div class="instructions">
            <h2>Tips & Tricks</h2>
            <ul>
                <li>Use the <strong>slider</strong> to navigate through MRI.
                Keyboard can also be used after clicking on slider.</li>
                <li>Check the <strong>Red</strong> crosshair which represents
                the placement of a given landmark.</li>
                <li>Press <strong>R</strong> to toggle between the subject
                scan and the protocol-defined placement (if one exists).</li>
                <li>Use the <strong>TAB</strong> key to navigate to
                next slider & MRI.</li>
            </ul>
        </div>
        <div class="viewer">
    """

    for label, images in subject_images.items():
        num_slices = len(images)
        has_reference = ( reference_images is not None
                         and label in reference_images )

        html_content += f"""
        <div class="container">
            <h2>Landmark: {afids_labels[label]}</h2>
            <img id="{label}_image" class="image"
            src="data:image/png;base64,{images[0]}"
            alt="MRI Slice">
            <input type="range" min="0"
            max="{num_slices - 1}" value="0"
            class="slider" id="{label}_slider">
        </div>
        """

    html_content += """
        </div>
        <script>
            const subjects = {};
    """

    for label, images in subject_images.items():
        has_reference = (
            reference_images is not None
            and label in reference_images
              )
        html_content += f"""
            subjects["{label}"] = {{
                targetImages: {images},
                referenceImages: {
                    reference_images[label] if has_reference else "null"
                    },
                showingReference: false
            }};
        """

    html_content += """
            document.addEventListener('DOMContentLoaded', () => {
                for (const [label, data] of Object.entries(subjects)) {
                    const slider = document.getElementById(`${label}_slider`);
                    const image = document.getElementById(`${label}_image`);

                    slider.addEventListener('input', () => updateImage(label));
                    document.addEventListener('keydown', (event) => {
                        if (event.key.toLowerCase() === 'r'
                            && data.referenceImages) {
                            data.showingReference = !data.showingReference;
                            updateImage(label);
                        }
                    });

                    function updateImage(label) {
                        const sliceIndex = document.getElementById(
                        `${label}_slider`
                        ).value;
                        const imageArray = subjects[label].showingReference
                            ? subjects[label].referenceImages
                            : subjects[label].targetImages;
                        document.getElementById(`${label}_image`).src =
                        'data:image/png;base64,' + imageArray[sliceIndex];
                    }
                }
            });
        </script>
    </body>
    </html>
    """

    with open(output_html, "w") as f:
        f.write(html_content)


def generate_interactive_mri_html(
    nii_path,
    fcsv_path,
    labels,
    ref_nii_path=None,
    ref_fcsv_path=None,
    jitter=2,
    zoom_radius=20,
    out_file_prefix="mri_viewer",
):
    """
    Generates an interactive HTML viewer for
    MRI slices based on fiducial coordinates for a single subject.

    Parameters:
    - nii_path: str
        Full path to the subject's NIfTI (.nii.gz) file.
    - fcsv_path: str
        Full path to the subject's FCSV (.fcsv) file.
    - labels: list of int
        Landmark indices to extract for the subject.
    - ref_nii_path: str or None
        Full path to reference NIfTI file (optional).
    - ref_fcsv_path: str or None
        Full path to reference FCSV file (optional).
    - jitter: int, optional (default=2)
        The number of pixels to expand around
        the coordinate in each slice for visualization.
    - zoom_radius: int, optional (default=20)
        The radius (in pixels) around the coordinate to extract for display.
    - out_file_prefix: str, optional (default="mri_viewer")
        The prefix for the output HTML file name.

    Notes:
    - Extracts specified fiducial coordinates from the subject's .fcsv file.
    - Maps the coordinates from world space to voxel
        space using the NIfTI affine transformation.
    - MRI slices centered around these coordinates are saved as images.
    - If a reference MRI and coordinate file are provided,
        the same process is applied.
    - Generates an interactive HTML viewer allowing keypress navigation.

    Returns:
    - None (outputs an HTML file with the
        interactive MRI viewer for the subject).
    """
    subject_key = os.path.basename(nii_path).replace(".nii.gz", "")
    target_img = nib.load(nii_path, mmap=True)
    affine_inv = np.linalg.inv(target_img.affine)
    target_data = target_img.get_fdata(dtype=np.float32)

    target_images = {}  # Dictionary to store images per landmark

    for label in labels:
        target_coord = extract_coordinates_from_fcsv(fcsv_path, label)
        if not target_coord:
            print(
                f"Coordinates for label '{label}'"
                + f"not found in subject '{subject_key}'."
            )
            continue

        target_voxel = np.round(affine_inv.dot((*target_coord, 1))).astype(int)
        target_images[label] = save_mri_slices_as_images(
            target_data, *target_voxel[:3], jitter, zoom_radius
        )

    reference_images = None
    if ref_nii_path and ref_fcsv_path:
        ref_img = nib.load(ref_nii_path, mmap=True)
        ref_data = ref_img.get_fdata(dtype=np.float32)
        reference_images = {}

        for label in labels:
            ref_coord = extract_coordinates_from_fcsv(ref_fcsv_path, label)
            if ref_coord:
                ref_voxel = np.round(
                    np.linalg.inv(ref_img.affine).dot((*ref_coord, 1))
                ).astype(int)
                reference_images[label] = save_mri_slices_as_images(
                    ref_data, *ref_voxel[:3], jitter, zoom_radius
                )

    out_file = f"{out_file_prefix}"

    # Pass the correct image dictionaries
    generate_html_with_keypress(target_images, reference_images, out_file)


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


def afids_to_fcsv(
    afid_coords: dict[int, NDArray],
    fcsv_output: os.PathLike[str] | str,
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
    with Path(fcsv_output).open(
        "w", encoding="utf-8", newline=""
        ) as out_fcsv_file:
        for line in header:
            out_fcsv_file.write(line)
        writer = csv.DictWriter(out_fcsv_file, fieldnames=AFIDS_FIELDNAMES)
        for row in fcsv:
            writer.writerow(row)


def fids_to_fcsv(fids, fcsv_template, fcsv_output):
    # Read in fcsv template
    with open(fcsv_template) as f:
        fcsv = [line.strip() for line in f]

    # Loop over fiducials
    for fid in range(1, fids.shape[0] + 1):
        # Update fcsv, skipping header
        line_idx = fid + 2
        centroid_idx = fid - 1
        fcsv[line_idx] = fcsv[line_idx].replace(
            f"afid{fid}_x", str(fids[centroid_idx][0])
        )
        fcsv[line_idx] = fcsv[line_idx].replace(
            f"afid{fid}_y", str(fids[centroid_idx][1])
        )
        fcsv[line_idx] = fcsv[line_idx].replace(
            f"afid{fid}_z", str(fids[centroid_idx][2])
        )

    # Write output fcsv
    with open(fcsv_output, "w") as f:
        f.write("\n".join(line for line in fcsv))


def get_download_dir():
    if "AUTOAFIDS_CACHE_DIR" in os.environ.keys():
        download_dir = os.environ["AUTOAFIDS_CACHE_DIR"]
    else:
        # create local download dir if it doesn't exist
        dirs = AppDirs("autoafids", "khanlab")
        download_dir = dirs.user_cache_dir
    return download_dir
