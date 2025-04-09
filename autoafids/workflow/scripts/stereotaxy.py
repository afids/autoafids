import pickle
import re
import warnings

import numpy as np
import pandas as pd

# Suppress specific warnings after all imports
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

# hardcoding right and left hemisphere afids for reflecting
right_afids = [
    "x_6",
    "x_8",
    "x_12",
    "x_15",
    "x_17",
    "x_21",
    "x_23",
    "x_25",
    "x_27",
    "x_29",
    "x_31",
    "y_6",
    "y_8",
    "y_12",
    "y_15",
    "y_17",
    "y_21",
    "y_23",
    "y_25",
    "y_27",
    "y_29",
    "y_31",
    "z_6",
    "z_8",
    "z_12",
    "z_15",
    "z_17",
    "z_21",
    "z_23",
    "z_25",
    "z_27",
    "z_29",
    "z_31",
]
left_afids = [
    "x_7",
    "x_9",
    "x_13",
    "x_16",
    "x_18",
    "x_22",
    "x_24",
    "x_26",
    "x_28",
    "x_30",
    "x_32",
    "y_7",
    "y_9",
    "y_13",
    "y_16",
    "y_18",
    "y_22",
    "y_24",
    "y_26",
    "y_28",
    "y_30",
    "y_32",
    "z_7",
    "z_9",
    "z_13",
    "z_16",
    "z_18",
    "z_22",
    "z_24",
    "z_26",
    "z_28",
    "z_30",
    "z_32",
]
combined_lables = [
    "AC",
    "PC",
    "ICS",
    "PMJ",
    "SIPF",
    "SLMS",
    "ILMS",
    "CUL",
    "IMS",
    "MB",
    "PG",
    "LVAC",
    "LVPC",
    "GENU",
    "SPLE",
    "ALTH",
    "SAMTH",
    "IAMTH",
    "IGO",
    "VOH",
    "OSF",
]
combined_lables = [
    element + axis for axis in ["x", "y", "z"] for element in combined_lables
]
right = [6, 8, 12, 15, 17, 21, 23, 25, 27, 29, 31]
left = [7, 9, 13, 16, 18, 22, 24, 26, 28, 30, 32]

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

def make_zero(num, threshold=0.0001):
    if abs(num) < threshold:
        return 0
    else:
        return num

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
    try:
        subject_id = re.search(r"(sub-\w+)", fcsv_path).group(1)
    except Exception as e:
       print("no subject id found")

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

def model_pred(
    in_fcsv: str,
    model: str,
    midpoint: str,
    slicer_tfm: str,
    template_fcsv: str,
    target_mcp: str,
    target_native: str,
):
    """
    Generate model predictions for fiducial points
    and transform coordinates to native space.

    Parameters
    ----------
        in_fcsv :: str
            Path to the input fiducial CSV file.
        model :: str
            Path to the trained model (pickle file).
        midpoint :: str
            Midpoint transformation matrix for fiducial alignment.
        slicer_tfm :: str
            ACPC transformation matrix from Slicer.
        template_fcsv :: str
            Template fiducial file for output format.
        target_mcp :: str
            Path to save MCP-transformed coordinates.
        target_native :: str
            Path to save native space coordinates.

    Returns
    -------
        None
    """
    # Transform input fiducial data using the specified transformation matrix
    fcsvdf_xfm = transform_afids(in_fcsv, slicer_tfm, midpoint)
    xfm_txt = fcsvdf_xfm[1]  # Transformation matrix in array form
    df_sub = dftodfml(fcsvdf_xfm[0])[0]

    # Compute MCP (midpoint of the collicular plate)
    # and center the fiducials on the MCP
    df_sub_mcp, mcp = mcp_origin(df_sub)

    # Reflect left hemisphere fiducials onto the right hemisphere.
    # This works because the data has already been ACPC-aligned
    # and MCP-centered.
    df_sub_mcp_l = df_sub_mcp.copy()
    df_sub_mcp_l.loc[
        :, df_sub_mcp_l.columns.str.contains("x")
    ] *= -1  # Flip 'x' coordinates to mirror

    # Drop left hemisphere fiducials from the original
    # and right hemisphere fiducials from the mirrored copy.
    # This retains the midline points with their original signs.
    df_sub_mcp = df_sub_mcp.drop(left_afids, axis=1)
    df_sub_mcp_l = df_sub_mcp_l.drop(right_afids, axis=1)

    # Standardize column names for concatenation
    df_sub_mcp.columns = combined_lables
    df_sub_mcp_l.columns = combined_lables

    # Combine the original and mirrored dataframes into a single dataset
    df_sub_mcp = pd.concat([df_sub_mcp, df_sub_mcp_l], ignore_index=True)

    # Replace near-zero values with exact zero
    # to avoid floating-point precision issues
    num_cols = df_sub_mcp.select_dtypes(include="number")
    cols_to_modify = (num_cols > -0.0001).all() & (num_cols < 0.0001).all()

    df_sub_mcp.loc[:, cols_to_modify] = (
        df_sub_mcp.loc[:, cols_to_modify]
        .applymap(make_zero)
    )


    # Load the trained model components from the pickle file
    with open(model, "rb") as file:
        objects_dict = pickle.load(file)

    # Extract preprocessing objects and Ridge regression models
    standard_scaler = objects_dict["standard_scaler"]
    pca = objects_dict["pca"]
    ridge_inference = [objects_dict["x"], objects_dict["y"], objects_dict["z"]]

    # Apply standard scaling and PCA transformation to the data
    df_sub_mcp = standard_scaler.transform(df_sub_mcp.values)
    df_sub_mcp = pca.transform(df_sub_mcp)

    # Make predictions using Ridge regression models for x, y, z coordinates
    y_sub = np.column_stack(
        [
            ridge.predict(df_sub_mcp) for ridge in ridge_inference
        ]
        )

    # Adjust the second predicted x-coordinate to reflect the left hemisphere
    y_sub[1, 0] *= -1

    # Save the predicted MCP-centered coordinates to a CSV file
    fids_to_fcsv(y_sub, template_fcsv, target_mcp)

    # Convert MCP-centered coordinates to native space
    stn_r_mcp = y_sub[0, :] + mcp.ravel()
    stn_l_mcp = y_sub[1, :] + mcp.ravel()

    # Create vectors for right and left fiducials with homogeneous coordinates
    vecr = np.hstack([stn_r_mcp.ravel(), 1])
    vecl = np.hstack([stn_l_mcp.ravel(), 1])

    # Apply the inverse transformation matrix
    # to convert coordinates to native space
    stn_r_native = np.linalg.inv(xfm_txt) @ vecr.T
    stn_l_native = np.linalg.inv(xfm_txt) @ vecl.T

    # Store the final native-space coordinates in a matrix
    stncoords = np.zeros((2, 3))
    stncoords[0, :] = stn_r_native[:3]
    stncoords[1, :] = stn_l_native[:3]

    # Save the native-space coordinates to the output file
    fids_to_fcsv(stncoords, template_fcsv, target_native)


model_pred(
    in_fcsv=snakemake.input.afidfcsv,
    model=snakemake.params.model,
    midpoint=snakemake.params.midpoint,
    slicer_tfm=snakemake.output.ACPC_txt,
    template_fcsv=snakemake.params.target_fcsv,
    target_mcp=snakemake.output.fcsv_mcp,
    target_native=snakemake.output.fcsv_native,
)
