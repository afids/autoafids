#!/usr/bin/env python
# coding: utf-8
import os  
import numpy as np  
import pandas as pd  
import pickle  
from pathlib import Path  

from sklearn.preprocessing import StandardScaler 
from sklearn.decomposition import PCA  
from sklearn.linear_model import Ridge

from .utils import transform_afids, fids_to_fcsv, mcp_origin, dftodfml, make_zero

#hardcoding right and left hemisphere afids
right_afids = ['x_6','x_8','x_12','x_15','x_17','x_21','x_23','x_25','x_27','x_29','x_31','y_6','y_8','y_12','y_15','y_17','y_21','y_23','y_25','y_27','y_29','y_31', 'z_6','z_8','z_12','z_15','z_17','z_21','z_23','z_25','z_27','z_29','z_31']
left_afids = ['x_7','x_9','x_13','x_16','x_18','x_22','x_24','x_26','x_28','x_30','x_32','y_7','y_9','y_13','y_16','y_18','y_22','y_24','y_26','y_28','y_30','y_32', 'z_7','z_9','z_13','z_16','z_18','z_22','z_24','z_26','z_28','z_30','z_32']
combined_lables = ['AC','PC', 'ICS', 'PMJ','SIPF','SLMS','ILMS','CUL','IMS','MB','PG','LVAC','LVPC','GENU','SPLE','ALTH','SAMTH','IAMTH','IGO','VOH','OSF']
combined_lables = [element + axis for axis in ['x', 'y', 'z'] for element in combined_lables]
r = [6,8,12,15,17,21,23,25,27,29,31]
l = [7,9,13,16,18,22,24,26,28,30,32]


def model_pred (in_fcsv, model, midpoint, slicer_tfm, template_fcsv, target_mcp, target_native):

    fcsvdf_xfm = transform_afids(in_fcsv, slicer_tfm, midpoint)
    xfm_txt = fcsvdf_xfm[1]
    df_sub = dftodfml(fcsvdf_xfm[0])
    df_sub_mcp , mcp = mcp_origin(df_sub) #compute MCP and center on MCP

    #reflect left side to the right (only works because coordiante data has been already ACPC transformed and MCP centered)
    df_sub_mcp_l = df_sub_mcp.copy()
    for column in df_sub_mcp_l.columns:
        if 'x' in column:
            df_sub_mcp_l[column] *= -1

    #drop left labels from right df vice versa; note now "midline" points have a sign depending on where they were placed (i.e., little to the right vs left). Hypothesis is that this information may be meaningful and should NOT be controlled
    df_sub_mcp = df_sub_mcp.drop(left_afids,axis =1 )
    df_sub_mcp_l = df_sub_mcp_l.drop(right_afids,axis=1)
    df_sub_mcp.columns = combined_lables
    df_sub_mcp_l.columns = combined_lables
    combined = [df_sub_mcp, df_sub_mcp_l]

    df_sub_mcp = pd.concat(combined, ignore_index=True)
    cols_to_modify = (df_sub_mcp.select_dtypes(include='number') > -0.0001).all() & (df_sub_mcp.select_dtypes(include='number') < 0.0001).all()
    df_sub_mcp.loc[:, cols_to_modify] = df_sub_mcp.loc[:, cols_to_modify].applymap(make_zero)

    # Load the bundled objects dictionary
    with open(model, 'rb') as file:
        objects_dict = pickle.load(file)

    # Extract the objects from the dictionary
    standard_scaler = objects_dict['standard_scaler']
    pca = objects_dict['pca']
    ridge_inference = [objects_dict['x'],objects_dict['y'],objects_dict['z']]


    df_sub_mcp = standard_scaler.transform(df_sub_mcp.values)
    df_sub_mcp = pca.transform(df_sub_mcp)


    y_sub = np.column_stack([ridge.predict(df_sub_mcp) for ridge in ridge_inference]) #make prediction on transformed data

    #relative to MCP in acpc space
    stncoords = np.zeros((2,3)) #init a empty matrix for native space coordaintes

    y_sub[1,0] = y_sub[1,0] * -1

    fids_to_fcsv(y_sub, template_fcsv, target_mcp) #save model outputs

    #in native space
    stn_r_mcp = y_sub[0, :] + mcp.ravel()
    stn_l_mcp = y_sub[1, :] + mcp.ravel()

    vecr = np.hstack([stn_r_mcp.ravel(), 1])
    vecl = np.hstack([stn_l_mcp.ravel(), 1])


    stn_r_native = np.linalg.inv(xfm_txt) @ vecr.T #applying tranformations to coordiantes
    stn_l_native = np.linalg.inv(xfm_txt) @ vecl.T #applying tranformations to coordiantes 


    stncoords[0,:] = stn_l_native[:3]
    stncoords[1,:] = stn_r_native[:3]

    fids_to_fcsv(stncoords, template_fcsv,target_native) #save model outputs



if __name__ == "__main__":
    model_pred(
        in_fcsv=snakemake.input["afidfcsv"],
        model=snakemake.params["model"],
        midpoint=snakemake.params["midpoint"],
        slicer_tfm=snakemake.output["ACPC_txt"],
        template_fcsv= snakemake.params["target_template"],
        target_mcp= snakemake.output["fcsv_mcp"],
        target_native= snakemake.output["fcsv_native"],
        
    )