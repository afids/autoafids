#!/usr/bin/env python
# coding: utf-8

import nibabel as nib


def normalize(method, input_im, output_im):
    """

    To normalize/zscore MRI image.

    Parameters
    ----------
        method:: str
            Method to be incorporated for normalizing. There are 2 choices, either minmax or zscore.

        input_im:: str
            Input image that needs to be normalized.

        output_im:: str
            Name of the modified image

    Returns
    -------
        None

    """
    nii = nib.load(input_im)
    nii_affine = nii.affine
    nii_data = nii.get_fdata()

    if method == "minmax":
        nii_data_normalized = (nii_data - nii_data.min()) / (
            nii_data.max() - nii_data.min()
        )

    elif method == "zscore":
        nii_data_normalized = (nii_data - nii_data.mean()) / (nii_data.std())

    nib.save(nib.Nifti1Image(nii_data_normalized, affine=nii_affine), output_im)


if __name__ == "__main__":
    normalize(
        method=snakemake.params["norm_method"],
        input_im=snakemake.input["im_raw"],
        output_im=snakemake.output["im_norm"],
    )
