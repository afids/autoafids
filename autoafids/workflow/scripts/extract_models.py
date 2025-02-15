import json
import tarfile
import tempfile
from os import PathLike
from pathlib import Path

def extract_afids_model(
    model_path: PathLike[str] | str,
    out_path: PathLike[str] | str,
) -> Path:
    """
    Extract afids model

    Parameters
    ----------
        model_path : PathLike
            Path to tarfile containing model weights
        
        out_path :: str
            Output path for predicted fiducails

    Returns
    -------
        None
    """
    with tarfile.open(model_path, "r:gz") as tar_file:
        for member in tar_file.getmembers():
            tar_file.extractall(
                path=out_path,
                members=[
                    candidate
                    for candidate in tar_file.getmembers()
                    if candidate.name.startswith(f"{member.name}/")
                ],
            )

extract_afids_model(snakemake.input.model_path, snakemake.output.model_dir)