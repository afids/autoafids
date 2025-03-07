import os  
import numpy as np  
from argparse import ArgumentParser
import pandas as pd  
import pickle  
from pathlib import Path  

from sklearn.preprocessing import StandardScaler 
from sklearn.decomposition import PCA  
from sklearn.linear_model import Ridge

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

from utils import transform_afids, fids_to_fcsv, mcp_origin, dftodfml, make_zero

#hardcoding right and left hemisphere afids for reflecting
right_afids = ['x_6','x_8','x_12','x_15','x_17','x_21','x_23','x_25','x_27','x_29','x_31','y_6','y_8','y_12','y_15','y_17','y_21','y_23','y_25','y_27','y_29','y_31', 'z_6','z_8','z_12','z_15','z_17','z_21','z_23','z_25','z_27','z_29','z_31']
left_afids = ['x_7','x_9','x_13','x_16','x_18','x_22','x_24','x_26','x_28','x_30','x_32','y_7','y_9','y_13','y_16','y_18','y_22','y_24','y_26','y_28','y_30','y_32', 'z_7','z_9','z_13','z_16','z_18','z_22','z_24','z_26','z_28','z_30','z_32']
combined_lables = ['AC','PC', 'ICS', 'PMJ','SIPF','SLMS','ILMS','CUL','IMS','MB','PG','LVAC','LVPC','GENU','SPLE','ALTH','SAMTH','IAMTH','IGO','VOH','OSF']
combined_lables = [element + axis for axis in ['x', 'y', 'z'] for element in combined_lables]
r = [6,8,12,15,17,21,23,25,27,29,31]
l = [7,9,13,16,18,22,24,26,28,30,32]

def model_pred(
    in_fcsv: str,
    model: str,
    midpoint: str,
    slicer_tfm: str,
    template_fcsv: str,
    target_mcp: str,
    target_native: str):

    """
    Generate model predictions for fiducial points and transform coordinates to native space.

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

    # Compute MCP (midpoint of the collicular plate) and center the fiducials on the MCP
    df_sub_mcp, mcp = mcp_origin(df_sub)

    # Reflect left hemisphere fiducials onto the right hemisphere.
    # This works because the data has already been ACPC-aligned and MCP-centered.
    df_sub_mcp_l = df_sub_mcp.copy()
    df_sub_mcp_l.loc[:, df_sub_mcp_l.columns.str.contains('x')] *= -1  # Flip 'x' coordinates to mirror

    # Drop left hemisphere fiducials from the original and right hemisphere fiducials from the mirrored copy.
    # This retains the midline points with their original signs.
    df_sub_mcp = df_sub_mcp.drop(left_afids, axis=1)
    df_sub_mcp_l = df_sub_mcp_l.drop(right_afids, axis=1)

    # Standardize column names for concatenation
    df_sub_mcp.columns = combined_lables
    df_sub_mcp_l.columns = combined_lables

    # Combine the original and mirrored dataframes into a single dataset
    df_sub_mcp = pd.concat([df_sub_mcp, df_sub_mcp_l], ignore_index=True)

    # Replace near-zero values with exact zero to avoid floating-point precision issues
    cols_to_modify = (
        (df_sub_mcp.select_dtypes(include='number') > -0.0001).all() & 
        (df_sub_mcp.select_dtypes(include='number') < 0.0001).all()
    )
    df_sub_mcp.loc[:, cols_to_modify] = df_sub_mcp.loc[:, cols_to_modify].applymap(make_zero)

    # Load the trained model components from the pickle file
    with open(model, 'rb') as file:
        objects_dict = pickle.load(file)

    # Extract preprocessing objects and Ridge regression models
    standard_scaler = objects_dict['standard_scaler']
    pca = objects_dict['pca']
    ridge_inference = [objects_dict['x'], objects_dict['y'], objects_dict['z']]

    # Apply standard scaling and PCA transformation to the data
    df_sub_mcp = standard_scaler.transform(df_sub_mcp.values)
    df_sub_mcp = pca.transform(df_sub_mcp)

    # Make predictions using Ridge regression models for x, y, z coordinates
    y_sub = np.column_stack([ridge.predict(df_sub_mcp) for ridge in ridge_inference])

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

    # Apply the inverse transformation matrix to convert coordinates to native space
    stn_r_native = np.linalg.inv(xfm_txt) @ vecr.T
    stn_l_native = np.linalg.inv(xfm_txt) @ vecl.T

    # Store the final native-space coordinates in a matrix
    stncoords = np.zeros((2, 3))
    stncoords[0, :] = stn_r_native[:3]
    stncoords[1, :] = stn_l_native[:3]
    
    # Save the native-space coordinates to the output file
    fids_to_fcsv(stncoords, template_fcsv, target_native)

if __name__ == "__main__":
    model_pred(
        in_fcsv= snakemake.input["afidfcsv"],
        model= snakemake.params["model"],
        midpoint= snakemake.params["midpoint"],
        slicer_tfm= snakemake.output["ACPC_txt"],
        template_fcsv= snakemake.params["target_fcsv"],
        target_mcp= snakemake.output["fcsv_mcp"],
        target_native= snakemake.output["fcsv_native"]  
    )