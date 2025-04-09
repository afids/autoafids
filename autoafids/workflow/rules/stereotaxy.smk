stereotaxy_target = config["stereotaxy"]


rule afidspred:
    input:
        afidfcsv=bids(
            root=root,
            datatype="afids-cnn",
            desc="afidscnn",
            suffix="afids.fcsv",
            **inputs[config["modality"]].wildcards
        ),
    output:
        fcsv_native=bids(
            root=root,
            datatype="stereotaxy",
            desc=stereotaxy_target,
            suffix="native.fcsv",
            **inputs[config["modality"]].wildcards
        ),
        fcsv_mcp=bids(
            root=root,
            datatype="stereotaxy",
            desc=stereotaxy_target,
            suffix="mcp.fcsv",
            **inputs[config["modality"]].wildcards
        ),
        ACPC_txt=bids(
            root=work,
            datatype="ACPCtransforms",
            desc="transform",
            suffix="ACPC.txt",
            **inputs[config["modality"]].wildcards
        ),
    params:
        model=str(Path(workflow.basedir).parent / config[stereotaxy_target]),
        midpoint="PMJ",
        target_fcsv=str(
            Path(workflow.basedir).parent / config["template_fcsv"]
        ),
    conda:
        "../envs/skimage.yaml"
    script:
        "../scripts/stereotaxy.py"
