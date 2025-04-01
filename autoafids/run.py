#!/usr/bin/env python3
from pathlib import Path
from autoafids.workflow.lib import utils as utils
import os

from snakebids import bidsapp, plugins

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.Version(distribution="autoafids"),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
    ]
)

# Set the conda prefix directory
conda_prefix = str(utils.get_download_dir()) +  "/" + "conda"

# Set the environment variable SNAKEMAKE_CONDA_PREFIX
os.environ["SNAKEMAKE_CONDA_PREFIX"] = str(conda_prefix)

def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
