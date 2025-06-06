# Command Line Interface (CLI)

## Autoafids CLI
The following can also be seen by entering `autoafids -h` into your terminal.

These are all the required and optional arguments Autoafids accepts in order to 
run flexibly on many different input data types and with many options. In most 
cases, only the required arguments are needed.

```{argparse}
---
filename: ../autoafids/run.py
func: get_parser
prog: autoafids
---
```

## Snakemake CLI
In addition to the above command line arguments, Snakemake arguments can also be
passed at the Autoafids command line.

The most critical of these is the `--cores / -c` and `--force-output` arguments,
which are **required** arguments for Autoafids.

The complete list of [Snakemake](https://snakemake.readthedocs.io/en/stable/) 
arguments are below, and most act to determine your environment and app
behaviours. They will likely only need to be used for running in cloud
environments or troubleshooting. These can be listed from the command line with
`autoafids --help-snakemake`.

```{argparse}
---
module: snakemake
func: get_argument_parser
prog: snakemake
---
```