# Output Files
In the specified `/path/to/output/dir`, there will be the following outputs:

```
/path/to/output/dir/
└── config
├── logs
├── sub-{subject}
├── work
├── .snakebids
└── .snakemake
```

The `config` folder, along with the hidden `.snakebids` and `.snakemake` folders
contain a record of the code and parameters used, and paths to the inputs.

## Work Directory 
After running the workflow, the `/path/to/output/dir` folder will contain a `work` directory. All the preprocessed nii.gz images for the autoafids CNN will be in the `work` directory with the following structure:

For a T1w image:
```
work/
└── sub-{subject}
    ├── n4biascorr
    ├── normalize
    ├── registration
    └── resample
```

For a T2w, FLAIR, or ct image:
```
work/
└── sub-{subject}
    ├── SynthSR
    ├── normalize
    ├── registration
    └── resample
```

When the --stereotaxy flag is specified, the ACPC tranforms are slso saved in the `work` directory:
```
work/
└── sub-{subject}
    └── ACPCtransforms
        └──*ACPC.txt
```

## sub-{subject} directory 
In the `sub-{subject}` directory an fcsv with the 32 anatomical fiducial locations is produced in the `afids-cnn` directory. When the --stereotaxy flag is used, an additional `stereotaxy` directory is produced with the stereotactic target locations for native space and mcp in fcsv files created. When the --fidqc flag is used, and additional `fidqc` directory is created with a `*afids.html` landmark qc file. 

```
sub-{subject}/
    ├── afids-cnn
    ├── fidqc
    └── stereotaxy
```
