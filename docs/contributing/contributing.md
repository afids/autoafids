# Contributing to Autoafids

Autoafids python package uses Poetry pacakge manager to manage its dependencies. You’ll need it installed on your machine before contributing to the software. Installation instructions can be found on the 
[Poetry website](https://python-poetry.org/docs/master/#installation).

Autoafids primarily caters to T1w modality images in which case it doesn't need to depend on pacakges outside of python but for other modalities it has a single dependency outside of python, i.e., `synthsr`. We strongly recommend using Autoafids with the `--use-singularity` flag in case of non T1w modalities, which will pull and use the required containers, unless you are comfortable installing and using `synthsr` tool yourself.

Note: These instructions are only recommended if you are making changes to the Autoafids codebase and committing these back to the repository or if you are using Snakemake’s cluster execution profiles.

## Setup the development environment

Once Poetry is available, clone this repository and install all dependencies (including `dev`):

```
git clone https://github.com/afids/autoafids.git
cd autoafids
poetry install --with dev 
```

Poetry will automatically create a virtual environment. To customize where 
these virtual environments are stored, see the poetry docs 
[here](https://python-poetry.org/docs/configuration/).

Then, you can run autoafids:

```
autoafids -h
```

or you can activate a virtual environment shell and run autoafids directly:

```
poetry shell
autoafids
```

You can exit the poetry shell with `exit`

## Running and fixing code format quality

Autoafids uses [poethepoet](https://github.com/nat-n/poethepoet) as a task runner.
You can see what commands are available by running:

```
poetry run poe 
```

We use a few tools, including `ruff`, `snakefmt`, and `yamlfix` to ensure 
formatting and style of our codebase is consistent. There are two task runners 
you can use to check and fix your code, which can be invoked with:

```
poetry run poe quality-check
poetry run poe quality
```

_Note: If you are in a poetry shell, you do not need to prepend `poetry run` to
the command._

## Dry-run / testing your workflow

Using Snakemake\'s dry-run option (`--dry-run`/`-n`) is an easy way to verify
any changes made to the workflow are working direcctly. The `tests/data` folder 
contains a _fake_ BIDS dataset (i.e. dataset with zero-sized files) that is 
useful for verifying different aspects of the workflow. These dry-run tests are 
part of the automated Github actions that are run for every commit.

You can invoke the pre-configured task via 
[poethepoet](https://github.com/nat-n/poethepoet) to perform a dry-run:

```
poetry run poe test_base
```

This performs a number of tests, involving different scenarios in which a user
may use autoafids.