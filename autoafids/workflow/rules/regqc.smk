from pathlib import Path
import glob

lead_dbs_dir = config.get("LEAD_DBS_DIR", False)
fmriprep_dir = config.get("FMRIPREP_DIR", False)


def get_warp_path(subject):
    if lead_dbs_dir:
        trans_dir = (
            Path(lead_dbs_dir)
            / f"sub-{subject}"
            / "normalization"
            / "transformations"
        )
        pattern = f"sub-{subject}*_from-MNI*_to-anchorNative_desc-ants.nii.gz"
        matches = list(trans_dir.glob(pattern))
        if matches:
            return str(matches[0])
        return None

    elif fmriprep_dir:
        trans_dir = Path(fmriprep_dir) / f"sub-{subject}" / "anat"
        pattern = f"sub-{subject}*from-T1w_to-MNI*_mode-image_xfm.h5"
        matches = list(trans_dir.glob(pattern))
        if matches:
            return str(matches[0])
        return None

    else:
        raise ValueError(
            "No LEAD-DBS or fMRIPrep directory provided for warp."
        )


def get_optional_matrix_path(subject):
    if lead_dbs_dir:
        trans_dir = (
            Path(lead_dbs_dir)
            / f"sub-{subject}"
            / "coregistration"
            / "transformations"
        )
        pattern = f"sub-{subject}_desc-precoreg_*T1w.mat"
        matches = list(trans_dir.glob(pattern))
        if matches:
            return str(matches[0])
        return []
    if fmriprep_dir:
        return []


def get_resampled_im(subject):
    if lead_dbs_dir:
        pattern = str(
            Path(lead_dbs_dir)
            / f"sub-{subject}"
            / "normalization"
            / "anat"
            / f"sub-{subject}_ses-preop_space-*_desc-preproc_acq-*_T1w.nii*"
        )
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
        else:
            raise FileNotFoundError(
                f"No resampled image found for subject {subject} in LEAD-DBS"
            )

    elif fmriprep_dir:
        pattern = str(
            Path(fmriprep_dir)
            / f"sub-{subject}"
            / "anat"
            / f"sub-{subject}*_space-MNI*_desc-preproc_T1w.nii.gz"
        )
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
        else:
            raise FileNotFoundError(
                f"No resampled image found for subject {subject} in fMRIPrep"
            )

    else:
        raise ValueError(
            "No LEAD-DBS or fMRIPrep directory provided for resampled image."
        )


def get_ref_paths():
    if lead_dbs_dir:
        refimage = str(
            Path(workflow.basedir).parent / config["templatet1w_lead"]
        )
        refcoordinate = str(
            Path(workflow.basedir).parent / config["fcsv_mni_lead"]
        )
    else:
        refimage = str(Path(workflow.basedir).parent / config["templatet1w"])
        refcoordinate = str(Path(workflow.basedir).parent / config["fcsv_mni"])
    return refimage, refcoordinate


rule regqc:
    input:
        afidfcsv=bids(
            root=root,
            datatype="afids-cnn",
            desc="afidscnn",
            suffix="afids.fcsv",
            **inputs[config["modality"]].wildcards
        ),
        im=lambda wildcards: get_resampled_im(wildcards.subject),
        warp=lambda wildcards: get_warp_path(wildcards.subject),
        optional_matrix=lambda wildcards: get_optional_matrix_path(
            wildcards.subject
        ),
    output:
        html=bids(
            root=root,
            datatype="regqc",
            desc="reg",
            suffix="qc.html",
            **inputs[config["modality"]].wildcards
        ),
        csv=bids(
            root=root,
            datatype="regqc",
            desc="reg",
            suffix="qc.csv",
            **inputs[config["modality"]].wildcards
        ),
        fcsv=bids(
            root=root,
            datatype="regqc",
            desc="reg",
            suffix="afids.fcsv",
            **inputs[config["modality"]].wildcards
        ),
    params:
        refim=lambda wildcards: get_ref_paths()[0],
        refcoord=lambda wildcards: get_ref_paths()[1],
    conda:
        "../envs/regqc.yaml"
    script:
        "../scripts/regqc.py"
