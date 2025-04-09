# Installation

## BIDS App for Automatic Anatomical Fiducials (AutoAFIDs)

### Requirements
- Docker (Windows/macOS/Linux) or Singularity (Linux)
- For those wishing to contribute or modify the code, `pip install` or `poetry install` are also available (Linux), but will still require Singularity for certain dependencies. See [Contributing to AutoAFIDs](https://github.com/afids/autoafids/blob/main/CONTRIBUTING.md).
- ⚠️ Warning: Local installation on Apple M1 is currently not supported. You can still run the pipeline by building a docker container, but it'll be an emulated amd64 container and fairly slow. 

### Notes
- Inputs to AutoAFIDs should be a BIDS-formatted dataset.
- The output includes fiducial markers that are automatically positioned at standardized anatomical landmarks.

## Docker on Windows/macOS/Linux

The AutoAFIDs BIDS App is available as a Docker container. Instructions can be found in the [Docker](https://github.com/afids/autoafids/tree/main/docker) documentation page.

### Pros
- Compatible with non-Linux systems, including Apple Silicon (M1, M2)
- All dependencies are contained within a single container

### Cons
- Typically not possible on shared machines
- Cannot use Snakemake cluster execution profiles
- Cannot edit code within the container

## Singularity Container

The same Docker container can also be used with Singularity (now Apptainer). Instructions can be found in the [Singularity](https://github.com/afids/autoafids/tree/main/docker) documentation page.

### Pros
- All dependencies are packaged within a single `.sif` file
- Compatible with shared systems where Singularity is installed

### Cons
- Cannot use Snakemake cluster execution profiles
- Cannot edit code within the container

## Python Environment with Singularity Dependencies

Instructions can be found in the [Contributing](https://github.com/afids/autoafids/blob/main/CONTRIBUTING.md) documentation page.

### Pros
- Provides flexibility to modify and develop the code

### Cons
- Requires a system with Singularity for handling external dependencies
