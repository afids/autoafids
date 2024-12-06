#!/usr/bin/env python
# coding: utf-8

from nilearn.image import resample_img
import numpy as np
import nibabel as nib


def resample_to_isotropic_resolution(
    nifti_image_path, isotropic_resolution, output_img
):
    """
    Resample a NIfTI MRI image to the specified isotropic resolution.

    Parameters
    ----------
        nifti_image_path:: str
            Path to the NIfTI MRI image file.

        isotropic_resolution:: float
            The desired isotropic resolution in mm.

        output_img :: str
            Name of the modified image

    Returns
    -------
        None

    """
    # Load the NIfTI image
    img = nib.load(nifti_image_path)

    # Get the original affine
    original_affine = np.eye(3, 3)

    # Create the target affine with the desired isotropic resolution
    target_affine = np.copy(original_affine)
    np.fill_diagonal(target_affine, isotropic_resolution / 100)

    # Resample the image
    resampled_image = resample_img(
        img, target_affine=target_affine, interpolation="continuous"
    )

    resampled_image.to_filename(output_img)


if __name__ == "__main__":
    resample_to_isotropic_resolution(
        nifti_image_path=snakemake.input["im_normed"],
        isotropic_resolution=snakemake.params["res"],
        output_img=snakemake.output["resam_im"],
    )
