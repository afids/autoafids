rule clone_synthsr:
    output:
        unzip_dir=directory(Path(download_dir) / "SynthSR"),
    log:
        bids(
            root="logs",
            suffix="clone_synthsr.log",
        ),
    params:
        download_dir=download_dir,
    script:
        "../scripts/clone_synthsr_repo.py"


rule synthsr:
    input:
        im=bids(
            root=str(Path(config["bids_dir"])),
            datatype="anat",
            suffix=f"{config['modality']}.nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        synthsr_dir=Path(download_dir) / "SynthSR",
    output:
        SynthSR=bids(
            root=work,
            datatype="SynthSR",
            desc=f"{config['modality']}-SynthSR",
            suffix="T1w.nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
    log:
        bids(
            root="logs",
            suffix="synthsr_output.txt",
            **inputs[config["modality"]].wildcards
        ),
    resources:
        mem_mb=10000,
    threads: max(1, workflow.cores // max(1, available_mem_mb // 10000))
    params:
        modality=config["modality"],
        download_dir=download_dir,
    conda:
        "../envs/synthsr.yaml"
    script:
        "../scripts/run_synthsr.py"
