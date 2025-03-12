# Contributing to Autoafids

Autoafids python package uses Poetry pacakge manager to manage its dependencies. You’ll need it installed on your machine before contributing to the software. Installation instructions can be found on the 
[Poetry website](https://python-poetry.org/docs/master/#installation).

Autoafids primarily caters to T1w modality images in which case it doesn't need to depend on pacakges outside of python but for other modalities it has a single dependency outside of python, i.e., `synthsr`. We strongly recommend using Autoafids with the `--use-singularity` flag in case of non T1w modalities, which will pull and use the required containers, unless you are comfortable installing and using `synthsr` tool yourself.

Note: These instructions are only recommended if you are making changes to the Autoafids codebase and committing these back to the repository or if you are using Snakemake’s cluster execution profiles.