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

_AFIDS = [f"{i:02d}" for i in range(1, 33)]  # ["01", "02", ..., "32"]

if config.get("enable_sequential_inference", False):
    rule applyfidmodel_all:
        """Run 5-patch inference for ALL 32 AFIDs using MNI-registered prior location"""
        input:
            t1w=lambda wildcards: (
                bids(
                    root=work,
                    datatype="normalize",
                    desc=chosen_norm_method,
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
                )
                if config["modality"] != "T1w"
                else bids(
                    root=work,
                    datatype="resample",
                    desc=chosen_norm_method,
                    res=config["res"],
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
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
        output:
            coords=[
                bids(
                    root=work,
                    datatype="afids-cnn",
                    afid=afid,
                    suffix="coord.txt",
                    **inputs[config["modality"]].wildcards,
                ) for afid in _AFIDS
            ],
            probs=[
                bids(
                    root=work,
                    datatype="afids-cnn",
                    afid=afid,
                    suffix="probmap.nii.gz",
                    **inputs[config["modality"]].wildcards,
                ) for afid in _AFIDS
            ],
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
                suffix="apply_all_prior.log",
                **inputs[config["modality"]].wildcards,
            ),
        params:
            ckpts=lambda wildcards: {
                key: str(Path(workflow.basedir).parent / path)
                for key, path in config["afids_inference"]["checkpoints"].items()
            }
        threads: 1

        conda:
            "../envs/pytorch.yaml"
        script:
            "../scripts/apply_with_prior_all.py"
else:
    # ── WITH PRIOR (SINGLE) ──────────────────────────────────────────────────────
    rule applyfidmodel_single:
        """Run 5-patch inference for ONE AFID using MNI-registered prior location."""
        input:
            t1w=lambda wildcards: (
                bids(
                    root=work,
                    datatype="normalize",
                    desc=chosen_norm_method,
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
                )
                if config["modality"] != "T1w"
                else bids(
                    root=work,
                    datatype="resample",
                    desc=chosen_norm_method,
                    res=config["res"],
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
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
        output:
            coord=bids(
                root=work,
                datatype="afids-cnn",
                afid="{afid}",
                suffix="coord.txt",
                **inputs[config["modality"]].wildcards,
            ),
            prob=bids(
                root=work,
                datatype="afids-cnn",
                afid="{afid}",
                suffix="probmap.nii.gz",
                **inputs[config["modality"]].wildcards,
            ),
        log:
            bids(
                root="logs",
                afid="{afid}",
                suffix="landmark.log",
                **inputs[config["modality"]].wildcards,
            ),
        wildcard_constraints:
            afid=r"\d{2}",
        params:
            ckpt_path=lambda wildcards: str(
                Path(workflow.basedir).parent
                / config["afids_inference"]["checkpoints"][f"afid_{int(wildcards.afid):02d}"]
            ),
        threads: 1

        conda:
            "../envs/pytorch.yaml"
        script:
            "../scripts/apply_with_prior_single.py"

    rule applyfidmodel_gather:
        """Collect all 32 per-AFID coord files and write combined FCSV."""
        input:
            coords=lambda wildcards: expand(
                bids(
                    root=work,
                    datatype="afids-cnn",
                    afid="{afid}",
                    suffix="coord.txt",
                    **{k: getattr(wildcards, k) for k in inputs[config["modality"]].wildcards},
                ),
                afid=_AFIDS,
            ),
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
                suffix="gather.log",
                **inputs[config["modality"]].wildcards,
            ),
        threads: 1
        conda:
            "../envs/pytorch.yaml"
        script:
            "../scripts/gather_afids.py"

# ---------------------------------------------------------------------------

if config.get("enable_sequential_inference", False):
    rule applyfidmodel_noprior_all:
        """Whole-volume sliding-window inference for ALL 32 AFIDs (no prior needed)"""
        input:
            t1w=lambda wildcards: (
                bids(
                    root=work,
                    datatype="normalize",
                    desc=chosen_norm_method,
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
                )
                if config["modality"] != "T1w"
                else bids(
                    root=work,
                    datatype="resample",
                    desc=chosen_norm_method,
                    res=config["res"],
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
                )
            ),
        output:
            coords=[
                bids(
                    root=work,
                    datatype="afids-cnn-noprior",
                    afid=afid,
                    suffix="coord.txt",
                    **inputs[config["modality"]].wildcards,
                ) for afid in _AFIDS
            ],
            probs=[
                bids(
                    root=work,
                    datatype="afids-cnn-noprior",
                    afid=afid,
                    suffix="probmap.nii.gz",
                    **inputs[config["modality"]].wildcards,
                ) for afid in _AFIDS
            ],
            fcsv=bids(
                root=root,
                datatype="afids-cnn",
                desc="afidscnn-noprior",
                suffix="afids.fcsv",
                **inputs[config["modality"]].wildcards,
            ),
        log:
            bids(
                root="logs",
                suffix="noprior-apply_all.log",
                **inputs[config["modality"]].wildcards,
            ),
        params:
            ckpts=lambda wildcards: {
                key: str(Path(workflow.basedir).parent / path)
                for key, path in config["afids_inference"]["checkpoints"].items()
            }
        threads: 1
        conda:
            "../envs/pytorch.yaml"
        script:
            "../scripts/apply_noprior_all.py"
else:
    # ── WITHOUT PRIOR (SINGLE) ──────────────────────────────────────────────────
    rule applyfidmodel_noprior_single:
        """Whole-volume sliding-window inference for ONE AFID — no prior needed."""
        input:
            t1w=lambda wildcards: (
                bids(
                    root=work,
                    datatype="normalize",
                    desc=chosen_norm_method,
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
                )
                if config["modality"] != "T1w"
                else bids(
                    root=work,
                    datatype="resample",
                    desc=chosen_norm_method,
                    res=config["res"],
                    suffix="T1w.nii.gz",
                    **inputs[config["modality"]].wildcards
                )
            ),
        output:
            coord=bids(
                root=work,
                datatype="afids-cnn-noprior",
                afid="{afid}",
                suffix="coord.txt",
                **inputs[config["modality"]].wildcards,
            ),
            prob=bids(
                root=work,
                datatype="afids-cnn-noprior",
                afid="{afid}",
                suffix="probmap.nii.gz",
                **inputs[config["modality"]].wildcards,
            ),
        log:
            bids(
                root="logs",
                afid="{afid}",
                suffix="noprior-landmark.log",
                **inputs[config["modality"]].wildcards,
            ),
        wildcard_constraints:
            afid=r"\d{2}",
        params:
            ckpt_path=lambda wildcards: str(
                Path(workflow.basedir).parent
                / config["afids_inference"]["checkpoints"][f"afid_{int(wildcards.afid):02d}"]
            ),
        threads: 1
        conda:
            "../envs/pytorch.yaml"
        script:
            "../scripts/apply_noprior_single.py"

    rule applyfidmodel_noprior_gather:
        """Collect all 32 no-prior coord files and write the combined FCSV."""
        input:
            coords=lambda wildcards: expand(
                bids(
                    root=work,
                    datatype="afids-cnn-noprior",
                    afid="{afid}",
                    suffix="coord.txt",
                    **{k: getattr(wildcards, k) for k in inputs[config["modality"]].wildcards},
                ),
                afid=_AFIDS,
            ),
        output:
            fcsv=bids(
                root=root,
                datatype="afids-cnn",
                desc="afidscnn-noprior",
                suffix="afids.fcsv",
                **inputs[config["modality"]].wildcards,
            ),
        log:
            bids(
                root="logs",
                suffix="noprior-gather.log",
                **inputs[config["modality"]].wildcards,
            ),
        threads: 1
        conda:
            "../envs/pytorch.yaml"
        script:
            "../scripts/gather_afids.py"