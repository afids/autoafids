include: "cnn.smk"

stereotaxy_target = config["stereotaxy"]

rule afidspred:
    input:
        afidfcsv=rules.gen_fcsv.output.fcsv
    output:
        fcsv_native=bids(
            root=root,
            datatype="stereotaxy",
            desc=stereotaxy_target,
            suffix="native.fcsv",
            **inputs["t1w"].wildcards
        ),
        fcsv_mcp=bids(
            root=root,
            datatype="stereotaxy",
            desc=stereotaxy_target,
            suffix="mcp.fcsv",
            **inputs["t1w"].wildcards
        ),        
        ACPC_txt=bids(
            root=work,
            datatype="ACPCtransforms",
            desc="transform",
            suffix="ACPC.txt",
            **inputs["t1w"].wildcards
        )
    params:
        model=str(Path(workflow.basedir).parent / config[stereotaxy_target]),
        midpoint="PMJ",
        target_fcsv = str(Path(workflow.basedir).parent / config['template_fcsv'])
    script:
        "../scripts/stereotaxy.py"
