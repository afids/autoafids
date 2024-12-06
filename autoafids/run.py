#!/usr/bin/env python3
from pathlib import Path

from snakebids import bidsapp, plugins

app = bidsapp.app(
    [
        plugins.SnakemakeBidsApp(Path(__file__).resolve().parent),
        plugins.Version(distribution="autoafids"),
    ]
)


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


if __name__ == "__main__":
    app.run()
