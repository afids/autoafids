#!/usr/bin/env python3
from pathlib import Path

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

if "__file__" not in globals():
    __file__ = "../autoafids/run.py"


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
