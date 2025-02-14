from __future__ import annotations

import numpy as np
import pandas as pd
import os
import re
from itertools import cycle
import random as rand
import base64
from io import BytesIO

import csv
from pathlib import Path

from numpy.typing import NDArray

# Dictionary for AFID labels
afids_labels = {
    1: 'AC', 2: 'PC', 3: 'ICS', 4: 'PMJ', 5: 'SIPF', 6: 'RSLMS',
    7: 'LSLMS', 8: 'RILMS', 9: 'LILMS', 10: 'CUL', 11: 'IMS', 12: 'RMB',
    13: 'LMB', 14: 'PG', 15: 'RLVAC', 16: 'LLVAC', 17: 'RLVPC', 18: 'LLVPC',
    19: 'GENU', 20: 'SPLE', 21: 'RALTH', 22: 'LALTH', 23: 'RSAMTH',
    24: 'LSAMTH', 25: 'RIAMTH', 26: 'LIAMTH', 27: 'RIGO', 28: 'LIGO',
    29: 'RVOH', 30: 'LVOH', 31: 'ROSF', 32: 'LOSF'
}


def fcsvtodf(fcsv_path):
    """
    Convert a .fcsv file (assumes RAS coordinate system) to a ML-friendly dataframe and return the cleaned xyz coordinates.
    
    Parameters:
    - fcsv_path: str, path to the .fcsv file
    
    Returns:
    - df_xyz_clean: pandas.DataFrame with cleaned x, y, z coordinates
    - num_points: int, number of fiducial points
    """

    # Extract the subject ID from the file path (naming is in bids-like)
    subject_id = re.search(r'(sub-\w+)', fcsv_path).group(1)

    # Read in .fcsv file, skip header
    df_raw = pd.read_table(fcsv_path,sep=',',header=2)

    #Extract the x, y, z coordiantes and store them in data science friendly format (i.e., features in cols and subject in rows)
    df_xyz = df_raw[['x','y','z']].melt().transpose()

    #Use number of row in fcsv to make number points
    colnames = [f'{axis}_{i % int(df_raw.shape[0]) + 1}' for axis in ['x', 'y', 'z'] for i in range(int(df_raw.shape[0]))]
    
    #Reassign features to be descriptive of coordinate
    df_xyz.columns = colnames

    #clean dataframe and pin correct subject name
    df_xyz_clean = df_xyz.drop('variable', axis= 0)
    df_xyz_clean = df_xyz_clean.rename(index={'value': subject_id})
    df_xyz_clean = df_xyz_clean.astype(float)

    return df_xyz_clean, df_raw.shape[0]

def dftodfml(fcsvdf):
    """
    Convert a datafrane (assumes RAS coordinate system) to a ML-friendly dataframe and return the cleaned xyz coordinates.
    
    Parameters:
    - fcsvdf: pandas.DataFrame
    
    Returns:
    - df_xyz_clean: pandas.DataFrame with cleaned x, y, z coordinates
    - num_points: int, number of fiducial points
    """

    #Extract the x, y, z coordiantes and store them in data science friendly format (i.e., features in cols and subject in rows)
    df_xyz = fcsvdf[['x','y','z']].melt().transpose()

    #Use number of row in fcsv to make number points
    colnames = [f'{axis}_{i % int(fcsvdf.shape[0]) + 1}' for axis in ['x', 'y', 'z'] for i in range(int(df_raw.shape[0]))]
    
    #Reassign features to be descriptive of coordinate
    df_xyz.columns = colnames

    #clean dataframe and pin correct subject name
    df_xyz_clean = df_xyz.drop('variable', axis= 0)
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
    coords1 = df[[f'x_{index1}', f'y_{index1}', f'z_{index1}']].to_numpy()
    coords2 = df[[f'x_{index2}', f'y_{index2}', f'z_{index2}']].to_numpy()

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
    coords1 = df[[f'x_{index1}', f'y_{index1}', f'z_{index1}']].to_numpy()
    coords2 = df[[f'x_{index2}', f'y_{index2}', f'z_{index2}']].to_numpy()

    # Compute the average as a numpy array
    xyz_average = (coords1 + coords2)/2

    return xyz_average.flatten()

def generate_slicer_file(matrix, filename):
    """
    Generate a .txt transformation file for 3D Slicer from a 4x4 matrix.

    Parameters: 
    - matrix: np.ndarray, 4x4 transformation matrix
    - output_path: str, path to store .txt file
    """
    D = np.array([
    [-1,  0,  0, 0],
    [ 0, -1,  0, 0],
    [ 0,  0,  1, 0],
    [ 0,  0,  0, 1]
    ])

    ras_inmatrix = np.linalg.inv(matrix)

    lps_inmatrix = D @ ras_inmatrix @ D

    # Extract rotation/scale and translation components
    rotation_scale = lps_inmatrix[0:3, 0:3].flatten()
    translation = lps_inmatrix[0:3, 3]

    # Format the content of the .tfm file
    tfm_content = "#Insight Transform File V1.0\n"
    tfm_content += "#Transform 0\n"
    tfm_content += "Transform: AffineTransform_double_3_3\n"
    tfm_content += "Parameters: " + " ".join(map(str, rotation_scale)) + " " + " ".join(map(str, translation)) + "\n"
    tfm_content += "FixedParameters: 0 0 0\n"

    # Write the content to the specified file
    with open(filename, 'w') as file:
        file.write(tfm_content)

def acpcmatrix(fcsv_path, midline, center_on_mcp = True, write_matrix = False, transform_file_name = None):
    """
    Computes a 4x4 transformation matrix aligning with the AC-PC axis.
    
    Parameters:
    - fcsv_path: str, path to the .fcsv file
    - midline: str or int, one of the 10 midline AFID points.
    - center_on_mcp: bool, If True, adds translation element to the ACPC matrix (centering on MCP).
    
    Returns:
    - matrix: np.ndarray, A 4x4 affine transformation matrix.
    """

    # A-P axis
    acpc = compute_distance(fcsv_path, 'AC', 'PC')  # vector from PC to AC
    yAxis = acpc[0] / acpc[1]  # unit vector defining anterior and posterior axis

    # R-L axis
    lataxis = compute_distance(fcsv_path, midline, 'AC')[0] # vector from AC to midline point 
    xAxis = np.cross(yAxis, lataxis) # vector defining left and right axis
    xAxis /= np.linalg.norm(xAxis) # unit vector defining left and right axis

    # S-I axis
    zAxis = np.cross(xAxis, yAxis)
    zAxis /= np.linalg.norm(zAxis)

    # Rotation matrix
    rotation = np.vstack([xAxis, yAxis, zAxis])
    translation = np.array([0, 0, 0])

    # Build 4x4 matrix
    matrix = np.eye(4)
    matrix[:3, :3] = rotation
    matrix[:3, 3] = translation

    # Orientation correction for matrix is midpoint placed below ACPC (i.e., at PMJ)
    if matrix[0][0] < 0:
        matrix = matrix * np.array([[-1,-1,-1,0],
                                    [1,1,1,0],
                                    [-1,-1,-1,0],
                                    [0,0,0,1]])

    if center_on_mcp: #need to compute MCP AFTER rotation is applied 
        mcp = compute_average(fcsv_path,'AC','PC') #MCP in native 
        matrix[:3, 3] = -np.dot(matrix[:3, :3], mcp) # MCP after rotation; negative because we set MCP to (0,0,0)

    if write_matrix: 
        generate_slicer_file(matrix, transform_file_name)

    return matrix

def mcp_origin(df_afids):
    """
    sets MCP as the origin point for each of the subjects
    """
    # extract MCP coordinates; defined as average point between AC and PC
    MCPx = (df_afids['x_1'] + df_afids['x_2'])/2
    MCPy = (df_afids['y_1'] + df_afids['y_2'])/2
    MCPz = (df_afids['z_1'] + df_afids['z_2'])/2

    # subtract MCP coordinates from afids at appropriate coords
    df_afids_mcpx = df_afids.transpose()[0:32] - MCPx
    df_afids_mcpy= df_afids.transpose()[32:64] - MCPy
    df_afids_mcpz= df_afids.transpose()[64:98] - MCPz

    #concat the three coords and take transpose
    frames = [df_afids_mcpx, df_afids_mcpy, df_afids_mcpz]
    df_afids_mcp = pd.concat(frames)
    df_afids_ori_mcp = df_afids_mcp.transpose()
    df_afids_ori_mcp = df_afids_ori_mcp.astype(float)
    
    return df_afids_ori_mcp , (np.array([MCPx.to_numpy(),MCPy.to_numpy(),MCPz.to_numpy()]))

def make_zero(num, threshold=0.0001):
    if abs(num) < threshold:
        return 0
    else:
        return num
    
def transform_afids(fcsv_path, slicer_tfm, midpoint):
    """
    Computes and applies an AC-PC transformation to AFID files and saves an fcsv with transfromed coordinates.

    Parameters:
    - fcsv_path: str, path to the .fcsv file for coordinates to be transformed
    - midline_fcsv: str, path to the .fcsv file which has the midline coordinates
    - output_dir: str, path of the directory to store transformed .fcsv file (if not specified, no output fcsv will be written)
    - midpoint: str or int, any midline AFID point
    - slicer_tfm: str, path to .txt file for slicer ACPC transform
   
    Returns:
    - tcoords: np.ndarray, transformed coordinates
    """
    
    # Compute the 4x4 AC-PC transformation matrix
    xfm_txt = acpcmatrix(fcsv_path, midpoint, slicer_tfm, center_on_mcp = False, write_matrix = True)

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
        raise ValueError("New coordinates do not match the number of rows in the original fcsv.")
    
    return fcsv_df,xfm_txt

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
    with Path(fcsv_output).open("w", encoding="utf-8", newline="") as out_fcsv_file:
        for line in header:
            out_fcsv_file.write(line)
        writer = csv.DictWriter(out_fcsv_file, fieldnames=AFIDS_FIELDNAMES)
        for row in fcsv:
            writer.writerow(row)

def fids_to_fcsv(fids, fcsv_template, fcsv_output):
    # Read in fcsv template
    with open(fcsv_template, "r") as f:
        fcsv = [line.strip() for line in f]

    # Loop over fiducials
    for fid in range(1, fids.shape[0]+1):
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