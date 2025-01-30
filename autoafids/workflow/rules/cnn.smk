# populate the AUTOAFIDS_CACHE_DIR folder as needed

def get_model():
    model_name = config["model"]

    local_model = config["resource_urls"].get(model_name, None)
    if local_model == None:
        print(f"ERROR: {model_name} does not exist.")

    return (Path(download_dir)/ "model" / Path(local_model).name).absolute()

rule download_cnn_model:
    params:
        url=config["resource_urls"][config["model"]],
        model_dir=Path(download_dir) / "model"
    output:
        model=get_model()
    shell:
        "mkdir -p {params.model_dir} && wget https://{params.url} -O {output}"

rule gen_fcsv:
    input:
        t1w = bids(
                root=str(Path(config["output_dir"]) / "resample"),
                datatype="anat",
                desc=chosen_norm_method,
                res=config["res"],
                suffix="T1w.nii.gz",
                **inputs["t1w"].wildcards,
                ), 
        model = get_model(),
        prior = bids(
            root=str(Path(config["output_dir"]) / "registration"),
            datatype="anat",
            space="native",
            desc="MNI",
            suffix="afids.fcsv",
            **inputs["t1w"].wildcards,
        ),
    output:
        fcsv = bids(
            root=str(Path(config["output_dir"]) / "afids-cnn"),
            desc="afidscnn",
            suffix="afids.fcsv",
            **inputs["t1w"].wildcards
        ),
    log:
        bids(
            root="logs",
            suffix="landmark.log",
            **inputs["t1w"].wildcards
        ),
    shell:
        'auto_afids_cnn_apply {input.t1w} {input.model} {output.fcsv} {input.prior}'