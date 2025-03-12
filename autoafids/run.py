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


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir"""
    app = bidsapp("../autoafids", skip_parse_args=True)
    add_dynamic_args(
        app.parser, app.config["parse_args"], app.config["pybids_inputs"]
    )
    return app.parser


if __name__ == "__main__":
    app.run()
