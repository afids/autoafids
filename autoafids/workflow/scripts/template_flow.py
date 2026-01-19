from templateflow import api as tflow

template_name = snakemake.params["template"]

if template_name in ['MNI305', 'MNIColin27']:
    template_path = tflow.get(
        str(template_name),
        desc=None,
        suffix='T1w',
        extension='nii.gz'
        )
else:
    template_path = tflow.get(
        str(template_name),
        desc=None,
        resolution=1,
        suffix='T1w',
        extension='nii.gz'
        )
