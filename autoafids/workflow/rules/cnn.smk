# populate the AUTOAFIDS_CACHE_DIR folder as needed


rule download_cnn_model:
    params:
        url=config["resource_urls"][config["model"]],
    output:
        unzip_dir=directory(Path(download_dir) / "models"),
    log:
        bids(
            root="logs",
            suffix="downlaod_model.log",
        ),
    shell:
        "wget https://{params.url} -O model.zip && "
        " unzip -q -d {output.unzip_dir} model.zip && "
        " rm model.zip"


rule applyfidmodel:
    input:
        t1w=lambda wildcards: (
            bids(
                root=work,
                datatype="normalize",
                desc=chosen_norm_method,
                suffix="T1w.nii.gz",
                **inputs[config["modality"]].wildcards,
            )
            if config["modality"] != "T1w"
            else bids(
                root=work,
                datatype="resample",
                desc=chosen_norm_method,
                res=config["res"],
                suffix="T1w.nii.gz",
                **inputs[config["modality"]].wildcards,
            )
        ),
        prior=bids(
            root=work,
            datatype="registration",
            space="native",
            desc="MNI",
            suffix="afids.fcsv",
            **inputs[config["modality"]].wildcards,
        ),
        model_dir=Path(download_dir) / "models",
    output:
        fcsv=bids(
            root=root,
            datatype="afids-cnn",
            desc="afidscnn",
            suffix="afids.fcsv",
            **inputs[config["modality"]].wildcards,
        ),
    log:
        bids(
            root="logs",
            suffix="landmark.log",
            **inputs[config["modality"]].wildcards,
        ),
    conda:
        "../envs/tensorflow.yaml"
    script:
        "../scripts/apply.py"
