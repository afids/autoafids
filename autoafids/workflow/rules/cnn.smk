# populate the AUTOAFIDS_CACHE_DIR folder as needed

def get_model():
    model_name = config["model"]

    local_model = config["resource_urls"].get(model_name, None)
    if local_model == None:
        print(f"ERROR: {model_name} does not exist.")

    return (Path(download_dir)/ "model" / Path(local_model).name).absolute()

def get_extracted_models():
    return (Path(download_dir)/ "extracted_model/").absolute()

rule download_cnn_model:
    params:
        url=config["resource_urls"][config["model"]],
        model_dir=Path(download_dir) / "model"
    output:
        model=get_model()
    shell:
        "mkdir -p {params.model_dir} && wget https://{params.url} -O {output}"

rule extract_models:
    input:
        model_path = get_model()
    output:
        model_dir = directory(get_extracted_models())
    script:
        "../scripts/extract_models.py"
    
rule gen_fcsv:
    input:
        t1w = bids(
                root=work,
                datatype="resample",
                desc=chosen_norm_method,
                res=config["res"],
                suffix="T1w.nii.gz",
                **inputs["t1w"].wildcards,
                ),
        extracted_model = get_extracted_models(),
        prior = bids(
            root=work,
            datatype="registration",
            space="native",
            desc="MNI",
            suffix="afids.fcsv",
            **inputs["t1w"].wildcards,
        ),
    params:
        model_path = get_model()
    output:
        fcsv = bids(
            root=root,
            datatype="afids-cnn",
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
        'auto_afids_cnn_apply {input.t1w} {params.model_path} {input.extracted_model} {output.fcsv} {input.prior}'