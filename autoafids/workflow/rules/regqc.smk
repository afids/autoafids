from pathlib import Path
import glob

lead_dbs_dir = config.get("LEAD_DBS_DIR", False)
fmriprep_dir = config.get("FMRIPREP_DIR", False)

template_name = config.get("template_flow", "")

template_dict = {
    "Agile12v2016": "resources/regqc/tpl-Agile12v2016_desc-groundtruth_afids.fcsv",
    "MNI152Lin": "resources/regqc/tpl-MNI152Lin_res-01_desc-groundtruth_afids.fcsv",
    "MNI152NLin2009cAsym": "resources/regqc/tpl-MNI152NLin2009cAsym_res-01_desc-groundtruth_afids.fcsv",
    "MNI152NLin6Sym": "resources/regqc/tpl-MNI152NLin6Asym_res-01_desc-groundtruth_afids.fcsv",
    "OASIS30ANTs": "resources/regqc/tpl-OASIS30ANTs_res-01_desc-groundtruth_afids.fcsv",
    "BigBrainSym": "resources/regqc/BigBrain_space-ICBM2009sym_desc-groundtruth_afids.fcsv",
    "MNI152NLin2009bAsym": "resources/regqc/tpl-MNI152NLin2009bAsym_res-1_desc-groundtruth_afids.fcsv",
    "MNI152NLin2009cSym": "resources/regqc/tpl-MNI152NLin2009cSym_res-1_desc-groundtruth_afids.fcsv",
    "MNI305": "resources/regqc/tpl-MNI305_desc-groundtruth_afids.fcsv",
    "PD25": "resources/regqc/tpl-PD25-res-01_desc-groundtruth_afids.fcsv",
    "fsaverage": "resources/regqc/tpl-fsaverage_res-01_den-41k_desc-groundtruth_afids.fcsv",
    "MNI152NLin2009bSym": "resources/regqc/tpl-MNI152NLin2009bSym_res-1_desc-groundtruth_afids.fcsv",
    "MNI152NLin6Asym": "resources/regqc/tpl-MNI152NLin6Asym_res-01_desc-groundtruth_afids.fcsv",
    "MNIColin27": "resources/regqc/tpl-MNIColin27_desc-groundtruth_afids.fcsv",
}


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
        pattern = f"sub-{subject}*from-MNI*_to-T1w_mode-image_xfm.h5"
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
    if template_name:
        refimage = directory(Path(download_dir) / "templateflow")

        refcoordinate = str(
            Path(workflow.basedir).parent / template_dict[template_name]
        )
    else:
        if lead_dbs_dir:
            refimage = str(
                Path(workflow.basedir).parent / config["templatet1w_lead"]
            )
            refcoordinate = str(
                Path(workflow.basedir).parent / config["fcsv_mni_lead"]
            )
        else:
            refimage = str(
                Path(workflow.basedir).parent / config["templatet1w"]
            )
            refcoordinate = str(
                Path(workflow.basedir).parent / config["fcsv_mni"]
            )
    return refimage, refcoordinate


rule download_template:
    params:
        template=template_name,
    output:
        template_path=directory(Path(download_dir) / "templateflow"),
    conda:
        "../envs/templateflow.yaml"
    script:
        "../scripts/template_flow.py"


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
        refim_dir=lambda wildcards: get_ref_paths()[0],
        refcoord=lambda wildcards: get_ref_paths()[1],
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
        template=template_name,
    conda:
        "../envs/regqc.yaml"
    script:
        "../scripts/regqc.py"
