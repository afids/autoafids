#!/usr/bin/env python3

from __future__ import annotations

# Forces TensorFlow to use CPU only 
import os
os.environ["CUDA_VISIBLE_DEVICES"] = ""

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
from tensorflow import keras

from .utils import afids_to_fcsv

logger = logging.getLogger(__name__)


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
        NDArray
            (fill description)
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
            (fill description)

        nii_affine :: NDArray
            (fill description)

    Returns
    -------
        fid_world :: NDArray
            (fill description)
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
            (fill description)

        nii_affine :: NDArray
            (fill description)
        
    Returns
    -------
        fid_voxel :: NDArray
            (fill description)
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
            (fill description)

        radius :: NDArray
            (fill description)

    Returns
    -------
        (fill Returns)
    """
    return tuple(slice(coord - radius, coord + radius + 1) for coord in centre[:3])


def slice_img(img: NDArray, centre: NDArray, radius: int) -> NDArray:
    """
    Slice images

    Parameters
    ----------
        img :: NDArray
            (fill description)
        
        centre :: NDArray
            (fill description)

        radius :: int
            (fill description)

    Returns
    -------
        NDArray 
            (fill desciption)

    """
    slices = gen_patch_slices(centre, radius)
    return img[slices[0], slices[1], slices[2]]


def predict_distances(
    radius: int,
    model: keras.model,
    mni_fid: NDArray,
    img: NDArray,
) -> NDArray:
    """
    Predict distances (elaborate)

    Parameters
    ----------
        radius :: int
            (fill description)

        model :: keras.model
            (fill description)

        mni_fid :: NDArray
            (fill description)
        
        img :: NDArray
            (fill description)

    Returns
    -------
        NDArray 
            (fill description)
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
            (fill description)

        img :: NDArray
            (fill description)

        mni_fid :: NDArray
            (fill description)

        radius :: int
            (fill description)

    Returns
    -------
        NDArray
            (fill description)
    """
    dim = (2 * radius) + 1
    print(f'min distance: {distances.min()}')
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
        logger.warning("No centroid found for this afid. Results may be suspect.")
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
    model: keras.model,
    radius: int,
    prior: PathLike[str] | str,
) -> NDArray:
    """
    Apply model

    Parameters
    ----------
        img :: nib.nifti1.Nifti1Image | nib.nifti1.Nifti1Pair
            (fill description)

        fid_label :: int
            (fill description)

        model :: keras.model
            (fill description)

        radius :: int
            (fill description)

        prior :: str 
            (fill description)

    Returns
    -------
        NDArray
            (fill description)
    """
    mni_fid_world = get_fid(load_fcsv(prior), fid_label)
    mni_img = img
    mni_fid_resampled = fid_world2voxel(
        mni_fid_world,
        mni_img.affine,
    )
    print('itr #1')
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
    print(f'itr #2')
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
    model_path: PathLike[str] | str,
    img: nib.nifti1.Nifti1Image | nib.nifti1.Nifti1Pair,
    prior: PathLike[str] | str,
) -> dict[int, NDArray]:
    """
    Apply all (elaborate)

    Parameters
    ----------
        model_path :: str
            (fill description)
        
        img :: nib.nifti1.Nifti1Image | nib.nifti1.Nifti1Pair
            (fill description)

        prior :: str
            (fill description)

    Returns
    -------
        afid_dict :: dict
            (fill description)
    """
    with tarfile.open(model_path, "r:gz") as tar_file:
        config_file = extract_config(tar_file)
        radius = int(json.load(config_file)["radius"])
        afid_dict: dict[int, NDArray] = {}
        for afid_label in range(1, 33):
            with tempfile.TemporaryDirectory() as model_dir:
                model = keras.models.load_model(
                    extract_afids_model(tar_file, model_dir, afid_label),
                )
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
            (fill description)

    Returns
    -------
        config_file :: IO[bytes]
            (fill description)
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
            (fill description)
        
        out_path :: str
            (fill description)

        afid_label :: int
            (fill description)

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
            (fill description)
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