# nnlm.smk  ─────────────────────────────────────────────────────────────────
# Rules for nnLandmark (nnLM) whole-volume AFID detection.
# Activated when `--detect_with_nnlm` is passed on the CLI.
#
# Rule graph:
#   download_nnlm_model  →  run_nnlm  →  nnlm_to_fcsv  →  rule all
# ─────────────────────────────────────────────────────────────────────────────

NNLM_MODEL_DIR = Path(download_dir) / "models" / "nnlm"


# ── 0. Download trained model ────────────────────────────────────────────────
rule download_nnlm_model:
    """Download the pre-trained nnLM model zip and extract it to the cache dir."""
    params:
        url=config["resource_urls"].get("nnlm", ""),
    output:
        model_dir=directory(NNLM_MODEL_DIR),
    log:
        bids(
            root="logs",
            suffix="download_nnlm_model.log",
        ),
    shell:
        "mkdir -p {output.model_dir} && "
        "wget -q {params.url} -O nnlm_model.zip > {log} 2>&1 && "
        "unzip -q -d {output.model_dir} nnlm_model.zip >> {log} 2>&1 && "
        "rm nnlm_model.zip"


# ── Helper: choose the right preprocessed T1w based on modality ──────────────
def _nnlm_t1w(wildcards):
    if config["modality"] != "T1w":
        return bids(
            root=work,
            datatype="normalize",
            desc=chosen_norm_method,
            suffix="T1w.nii.gz",
            **inputs[config["modality"]].wildcards,
        )
    return bids(
        root=work,
        datatype="resample",
        desc=chosen_norm_method,
        res=config["res"],
        suffix="T1w.nii.gz",
        **inputs[config["modality"]].wildcards,
    )


# ── 1. nnLM inference (all 32 AFIDs in one forward pass) ─────────────────────
rule run_nnlm:
    """Run nnLandmark whole-volume inference — all 32 AFIDs in a single forward pass."""
    input:
        t1w=_nnlm_t1w,
        model_dir=NNLM_MODEL_DIR,
    output:
        seg=bids(
            root=work,
            datatype="afids-nnlm",
            suffix="dseg.nii.gz",
            **inputs[config["modality"]].wildcards,
        ),
        coords_json=bids(
            root=work,
            datatype="afids-nnlm",
            suffix="coords.json",
            **inputs[config["modality"]].wildcards,
        ),
    params:
        dataset_id="800",
        fold=config.get("nnlm_fold", "0"),
        plans=config.get("nnlm_plans", "nnUNetResEncUNetMPlans"),
        checkpoint=config.get("nnlm_checkpoint", "checkpoint_final.pth"),
        device=config.get("nnlm_device", "cuda"),
        tmpdir=lambda wildcards: f"/tmp/nnlm_{wildcards.subject}",
    conda:
        "../envs/nnlm.yaml"
    threads: 4
    resources:
        gpus=lambda wildcards: 1 if config.get("nnlm_device", "cuda") == "cuda" else 0,
        mem_mb=16000,
    log:
        bids(
            root="logs",
            suffix="nnlm.log",
            **inputs[config["modality"]].wildcards,
        ),
    shell:
        """
        set -euo pipefail

        # Ensure log and output directories exist
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.seg})

        TMPDIR={params.tmpdir}
        mkdir -p $TMPDIR/input_dir $TMPDIR/output_dir

        # nnLM input naming convention: {{case}}_0000.nii.gz
        cp {input.t1w} $TMPDIR/input_dir/{wildcards.subject}_0000.nii.gz

        # Point nnLM env vars to the downloaded model
        export nnLM_raw=dummy
        export nnLM_preprocessed=dummy
        export nnLM_results={input.model_dir}

        python -m nnlandmark.inference.nnLandmark.predict_from_raw_data \
            -i $TMPDIR/input_dir \
            -o $TMPDIR/output_dir \
            -d {params.dataset_id} \
            -c 3d_fullres \
            -tr nnLandmark \
            -p {params.plans} \
            -f {params.fold} \
            -chk {params.checkpoint} \
            -device {params.device} \
            -npp 0 -nps 0 \
            > {log} 2>&1

        # Collect outputs
        cp $TMPDIR/output_dir/{wildcards.subject}.nii.gz {output.seg}
        cp $TMPDIR/output_dir/{wildcards.subject}.json    {output.coords_json}

        # Cleanup temp dir
        rm -rf $TMPDIR
        """


# ── 2. Convert JSON voxel coords → RAS world coords → Slicer FCSV ────────────
rule nnlm_to_fcsv:
    """Convert nnLM per-subject JSON (voxel space) to a Slicer-compatible FCSV (RAS mm)."""
    input:
        coords_json=bids(
            root=work,
            datatype="afids-nnlm",
            suffix="coords.json",
            **inputs[config["modality"]].wildcards,
        ),
        t1w=_nnlm_t1w,
    output:
        fcsv=bids(
            root=root,
            datatype="afids-cnn",
            desc="afidscnn-nnlm",
            suffix="afids.fcsv",
            **inputs[config["modality"]].wildcards,
        ),
    params:
        fcsv_template=str(Path(workflow.basedir).parent / "resources" / "dummy.fcsv")
    conda:
        "../envs/nibabel.yaml"
    script:
        "../scripts/nnlm_to_fcsv.py"
