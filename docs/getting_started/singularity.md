# Running Autoafids with Singularity

## Pre-requisites
1. Singularity / Apptainer is installed on your system. For more info, see
the detailed [Apptainer install instructions](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages).
1. The following command-line tools are installed:
    * wget
1. Sufficient disk-space (rough estimate)
    * in your `/tmp` folder (>30GB) to build the container
    * in your working folder to store the container (~20GB)
    * for AutoAFIDs outputs (~10GB per subject using default parameters)
1. Sufficient CPU and memory - the more you have, the faster it will run. We recommend at least 4 CPU cores 
   and 16GB memory if using default parameters.

## First time setup
Pull the container. This can be done from DockerHub, but requires a large 
amount of disk space in your `/tmp` folder, since it has to convert from a 
Docker container to a Singularity/Apptainer container. The example below pulls
the latest versioned container (replace `latest` with `vX.X.X` for a specific
version).

```
singularity pull docker://khanlab/autoafids:latest
```

_Note: If you encounter any errors pulling the container from DockerHub, it may
be because you are running out of disk space in your cache folders. You can 
change these locations by setting environment variables:_

```
export SINGULARITY_CACHEDIR=/YOURDIR/.cache/singularity
```

Run Autoafids without any arguments to print the short help:

```
singularity run -e autoafids_latest.sif
```

Use the `-h` option to get a detailed help listing:

```
singularity run -e autoafids_latest.sif -h
```

## Running an example

We will use the `tests` folder found from the 
[Github repository](https://github.com/afids/autoafids/tree/main/test/) to
demonstrate an example of how to run Autoafids:

```
singularity run -e autoafids_latest.sif tests/data/ tests/data/derivatives participant --force-output -n
```

### Explanation

The first three arguments to 
Autoafids (as with any BIDS App) are the input folder (`tests/data/`), the 
output folder (`tests/data/derivatives`), and the analysis level (`participant`). The `--force-output` 
flag enables writing output files even if the folder exists. The `--dry-run/-n` option lists what would run 
without actually running the workflow.

To see the shell commands for each rule:

```
singularity run -e autoafids_latest.sif tests/data/ tests/data/derivatives participant --force-output -np
```

To run the workflow with all available CPU cores:

```
singularity run -e autoafids_latest.sif tests/data/ tests/data/derivatives participant --force-output -p --cores all
```

_Note: You may need to bind paths to ensure proper file access:_

```
export SINGULARITY_BINDPATH=/data:/data
```

After completion, you will find results in your output folder, `tests/data/derivatives`.