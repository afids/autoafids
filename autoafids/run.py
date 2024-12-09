#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins

import logging

# Set up the root logger for general debug-level logging
# logging.basicConfig(level=logging.DEBUG)

# Prevent duplicate terminal log outputs from snakemake
snakemake_logger = logging.getLogger("snakemake.logging")
snakemake_logger.propagate = False  # Prevent duplicate logs
snakemake_logger.setLevel(logging.INFO)  # Adjust level if needed

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.BidsValidator(),
        plugins.Version(distribution="autoafids"),
        plugins.CliConfig("parse_args"),
        plugins.ComponentEdit("pybids_inputs"),
    ]
)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
