#!/usr/bin/env python3

from __future__ import annotations

# Forces TensorFlow to use CPU only 
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "" #only use CPU, makes generealization to any machine easier.
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"  # Suppresses INFO, WARNING, and ERROR logs


import json
import logging
import tarfile
import tempfile
from argparse import ArgumentParser
from os import PathLike
from pathlib import Path
from typing import IO
import nibabel as nib
import numpy as np
import pandas as pd
import skimage.measure
from numpy.typing import NDArray
import tensorflow as tf
from .utils import afids_to_fcsv


def load_fcsv(fcsv_path: PathLike[str] | str) -> pd.DataFrame:
    """
    Loads fcsv.

    Parameters
    ----------
        fcsv_path :: str
            Path to a fcsv file

    Returns
    -------
        pd.DataFrame
            A DataFrame containing the data from the FCSV file
    """
    return pd.read_csv(fcsv_path, sep=",", header=2)


# utils to factor out
def get_fid(fcsv_df: pd.DataFrame, fid_label: int) -> NDArray:
    """
    Extract specific fiducial's spatial coordinates.

    Parameters
    ----------
        fcsv_df :: pd.DataFrame
            Dataframe with the FCSV data.
        fid_label :: int
            Label (1-32) of the fiducial to grab.

    Returns
    -------
        NDArray :: np.array
            Array with x,y,z coordinates of selected fiducial
    """
    return fcsv_df.loc[fid_label - 1, ["x", "y", "z"]].to_numpy(
        dtype="single",
        copy=True,
    )


def fid_voxel2world(fid_voxel: NDArray, nii_affine: NDArray) -> NDArray:
    """
    Transform fiducials in voxel coordinates to world coordinates.
    
    Parameters
    ----------
        fid_voxel :: NDArray
            Voxel coordiantes for selected fiducial 

        nii_affine :: NDArray
            Image affine matrix for mapping between voxel to world coords 

    Returns
    -------
        fid_world :: NDArray
            Coordinates in world space 
    """
    translation = nii_affine[:3, 3]
    rotation = nii_affine[:3, :3]
    fid_world = rotation.dot(fid_voxel)+translation
    return fid_world.astype(float)

def fid_world2voxel(
    fid_world: NDArray,
    nii_affine: NDArray,
) -> NDArray:
    """
    Transform fiducials in world coordinates to voxel coordinates.
    
    Parameters
    ----------
        fid_world :: NDArray
            World coordiantes for selected fiducial

        nii_affine :: NDArray
            Image affine matrix for mapping between voxel to world coords 
        
    Returns
    -------
        fid_voxel :: NDArray
            Coordinates in voxel space 
    """
    inv_affine =  np.linalg.inv(nii_affine)
    translation = inv_affine[:3, 3]
    rotation = inv_affine[:3, :3]
    fid_voxel = rotation.dot(fid_world) + translation
    fid_voxel = np.rint(fid_voxel)
    return fid_voxel.astype(int)


def gen_patch_slices(centre: NDArray, radius: int) -> tuple[slice, slice, slice]:
    """
    Generate patch slices

    Parameters
    ----------
        centre :: NDArray
            Center of patch

        radius :: NDArray
            Radius around center of patch (model defined)

    Returns
    -------
        Image patches around center coordinate
    """
    return tuple(slice(coord - radius, coord + radius + 1) for coord in centre[:3])


def slice_img(img: NDArray, centre: NDArray, radius: int) -> NDArray:
    """
    Slice images

    Parameters
    ----------
        img :: NDArray
            input image to be predicted
        
        centre :: NDArray
            Center of patch

        radius :: int
            patch radius

    Returns
    -------
        NDArray 
            (fill desciption)

    """
    slices = gen_patch_slices(centre, radius)
    return img[slices[0], slices[1], slices[2]]


def predict_distances(
    radius: int,
    model: tf.keras.model,
    mni_fid: NDArray,
    img: NDArray,
) -> NDArray:
    """
    Predict distances (elaborate)

    Parameters
    ----------
        radius :: int
            patch radius

        model :: tf.keras.model
            ML model for predicting within patch 

        mni_fid :: NDArray
            MNI fid prior to limit scope of prediction
        
        img :: NDArray
            full image

    Returns
    -------
        NDArray 
            predicted label map 
    """
    dim = (2 * radius) + 1
    pred = np.reshape(slice_img(img, mni_fid, radius), (1, dim, dim, dim, 1))
    return model.predict(pred)


def process_distances(
    distances: NDArray,
    img: NDArray,
    mni_fid: NDArray,
    radius: int,
) -> NDArray:
    """
    Process distances (elaborate)

    Parameters
    ----------
        distances :: NDArray
            predicted distance map (i.e., ML model output)

        img :: NDArray
            input image to be predicted

        mni_fid :: NDArray
            MNI fid prior

        radius :: int
            Radius of configured model 

    Returns
    -------
        NDArray
            processed prediction to x,y,z coordinates
    """
    dim = (2 * radius) + 1
    arr_dis = np.reshape(distances[0], (dim, dim, dim))
    new_pred = np.full((img.shape), 100, dtype=float)
    slices = gen_patch_slices(mni_fid, radius)
    new_pred[slices[0], slices[1], slices[2]] = arr_dis
    transformed = np.exp(-0.5 * new_pred)
    thresh = np.percentile(transformed, 99)
    thresholded = transformed
    thresholded[thresholded < thresh] = 0
    thresholded = (thresholded * 1000000).astype(int)
    new = skimage.measure.regionprops(thresholded)
    if not new:
        print("No centroid found for this afid. Results may be suspect.")
        return np.array(
            np.unravel_index(
                np.argmax(transformed, axis=None),
                transformed.shape,
            ),
        )
    centroids = {
        key: [region.centroid[idx] for region in new]
        for idx, key in enumerate(["x", "y", "z"])
    }
    return np.array(
        [sum(centroids[key]) / len(centroids[key]) for key in ["x", "y", "z"]],
    )


def apply_model(
    img: nib.nifti1.Nifti1Image | nib.nifti1.Nifti1Pair,
    fid_label: int,
    model: tf.keras.model,
    radius: int,
    prior: PathLike[str] | str,
) -> NDArray:
    """
    Apply model

    Parameters
    ----------
        img :: nib.nifti1.Nifti1Image | nib.nifti1.Nifti1Pair
            input image to be predicted

        fid_label :: int
            fiducial label of interest as defined by the protocol 

        model :: tf.keras.model
            ML model for predicting within patch 

        radius :: int
            patch radius

        prior :: str 
            MNI fid prior fiducal file

    Returns
    -------
        NDArray
            processed prediction to x,y,z coordinates
    """
    mni_fid_world = get_fid(load_fcsv(prior), fid_label)
    mni_img = img
    mni_fid_resampled = fid_world2voxel(
        mni_fid_world,
        mni_img.affine,
    )
    img_data = img.get_fdata()
    distances = predict_distances(
        radius,
        model,
        mni_fid_resampled,
        img_data,
    )
    fid_resampled = process_distances(
        distances,
        img_data,
        mni_fid_resampled,
        radius,
    )
    #do it again to improve prediction
    fid_pred = np.rint(fid_resampled).astype(int)
    distances2 = predict_distances(
        radius,
        model,
        fid_pred,
        img_data,
    )
    fid_resampled2 = process_distances(
        distances2,
        img_data,
        fid_pred,
        radius,
    )
    return fid_voxel2world(fid_resampled2, img.affine)

def apply_all(
    model_path: PathLike,
    img: nib.nifti1.Nifti1Image,
    prior: PathLike,
) -> Dict[int, NDArray]:
    """
    Apply all models using a single architecture with dynamically loaded weights.

    Parameters
    ----------
        model_path : PathLike
            Path to tarfile containing model weights
        
        img : nib.nifti1.Nifti1Image
            Input image to be predicted

        prior : PathLike
            Path to MNI fid prior fiducial file

    Returns
    -------
        afid_dict : dict
            Dictionary of all predicted landmarks
    """
    with tarfile.open(model_path, "r:gz") as tar_file:
        # Extract configuration
        config_file = extract_config(tar_file)
        radius = int(json.load(config_file)["radius"])
        
        # Extract the shared model architecture
        with tempfile.TemporaryDirectory() as model_dir:
            model_architecture_path = extract_afids_model(tar_file, model_dir, 1)
            model = tf.keras.models.load_model(model_architecture_path)

        # Create an empty dictionary for storing predictions
        afid_dict: Dict[int, NDArray] = {}

        # Iterate through labels
        for afid_label in range(1, 33):
            print(f"Processing AFID Label: {afid_label}")

            # Locate the weight file for the current label directly in the tarfile
            weight_file_name = extract_afids_model(tar_file, model_dir, afid_label)
            
            
            # Load weights into memory (avoid disk I/O)
           # weight_file_data = BytesIO(tar_file.extractfile(weight_file_name).read())
            model.load_weights(weight_file_name).expect_partial()

            # Apply the model to make predictions
            afid_dict[afid_label] = apply_model(
                img,
                afid_label,
                model,
                radius,
                prior,
            )

    return afid_dict


def extract_config(tar_file: tarfile.TarFile) -> IO[bytes]:
    """
    Extract config (elaborate)

    Parameters
    ----------
        tar_file :: arfile.TarFile
            ML model .tar file genereated from train workflow 

    Returns
    -------
        config_file :: IO[bytes]
            ML model configration file 
    """
    try:
        config_file = tar_file.extractfile("config.json")
    except KeyError as err:
        missing_data = "config file"
        raise ArchiveMissingDataError(missing_data, tar_file) from err
    if not config_file:
        missing_data = "config file as file"
        raise ArchiveMissingDataError(missing_data, tar_file)
    return config_file


def extract_afids_model(
    tar_file: tarfile.TarFile,
    out_path: PathLike[str] | str,
    afid_label: int,
) -> Path:
    """
    Extract afids model

    Parameters
    ----------
        tar_file :: tarfile.TarFile
            ML model .tar file genereated from train workflow 
        
        out_path :: str
            Output path for predicted fiducails

        afid_label :: int
            fiducial label of interest as defined by the protocol 

    Returns
    -------
        None
    """
    for member in tar_file.getmembers():
        if member.isdir() and f"afid-{afid_label:02}" in member.name:
            tar_file.extractall(
                path=out_path,
                members=[
                    candidate
                    for candidate in tar_file.getmembers()
                    if candidate.name.startswith(f"{member.name}/")
                ],
            )

            return Path(out_path) / member.name
    msg = f"AFID {afid_label:02} model"
    raise ArchiveMissingDataError(msg, tar_file)


class ArchiveMissingDataError(Exception):
    def __init__(self, missing_data: str, tar_file: tarfile.TarFile) -> None:
        super().__init__(
            f"Required data {missing_data} not found in archive {tar_file}.",
        )


def gen_parser() -> ArgumentParser:
    """
    Gen parser (elaborate)

    Parameters
    ----------
        None
    
    Returns
    -------
        parser :: ArgumentParser
            parser to help navigate paths and configure functions
    """
    parser = ArgumentParser()
    parser.add_argument("img", help="The image for which to produce an FCSV.")
    parser.add_argument("model", help="The afids-CNN model to apply.")
    parser.add_argument("fcsv_path", help="The path to write the output FCSV.")
    parser.add_argument("sub_prior", help="The coordinates to define model prediction space")
    return parser


def main():
    """
    Main

    Parameters
    ----------
        None

    Returns
    -------
        None
    """
    args = gen_parser().parse_args()
    img = nib.nifti1.load(args.img)

    predictions = apply_all(args.model, img,args.sub_prior)
    afids_to_fcsv(predictions, args.fcsv_path)


if __name__ == "__main__":
    main()
