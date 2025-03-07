rule fidqc:
    input:
        afidfcsv = bids(
            root=root,
            datatype="afids-cnn",
            desc="afidscnn",
            suffix="afids.fcsv",
            **inputs[config['modality']].wildcards
        ),
        im=bids(
            root=str(Path(config["bids_dir"])),
            datatype="anat",
            suffix="T1w.nii.gz",
            **inputs[config['modality']].wildcards
        )if config['modality'] == "T1w" else
        bids(
            root=work,
            datatype="non_t1w_preprocessed",
            desc="non_t1w_preprocessed",
            suffix="T1w.nii.gz",
            **inputs[config['modality']].wildcards
        )
    output:
        html=bids(
            root=root,
            datatype="fidqc",
            desc="fidqc",
            suffix="afids.html",
            **inputs[config['modality']].wildcards
        )
    params:
        refim=str(Path(workflow.basedir).parent / config["templatet1w"]),
        refcoord=str(Path(workflow.basedir).parent / config["fcsv_mni"])
    script:
        '../scripts/qc.py'
