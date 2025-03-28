# Workflow Details

This section describes the Autoafids workflow (i.e. steps taken to produce 
intermediate and final files). Autoafids is a 
[Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow, and thus a 
directed acyclic graph (DAG) that is automatically configured based on a set of 
rules.

## Overall workflow

Below is an image exhibiting the workflow DAG. Each rounded rectangle underneath the modality input row (T1w, T2w, FLAIR, ct) in 
the DAG represents a rule (i.e. some code or script that produces an output), 
with arrows representing the connections (i.e. inputs / outputs) to these rules.

<img src="https://raw.githubusercontent.com/afids/autoafids/refs/heads/main/docs/dag.svg" width="800px">
