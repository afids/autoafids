import os
from templateflow import api as tflow 

template_name = snakemake.params["template"]
output_dir = snakemake.output["template_dir"]

os.environ['TEMPLATEFLOW_HOME'] = output_dir

tflow.get(template_name, desc=None, resolution=1, suffix='T1w', extension='nii.gz')