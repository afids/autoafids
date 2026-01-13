from templateflow import api as tflow

template_name = snakemake.params["template"]

template_path = tflow.get(str(template_name), desc=None, resolution=1, suffix='T1w', extension='nii.gz')
