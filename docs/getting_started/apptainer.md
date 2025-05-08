# Running AutoAFIDs with Apptainer (Singularity)

## Pre-requisities:
 1. Apptainer (or Singularity) is installed on your system. For more info, see the detailed [apptainer install instructions](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages).
 2. The following command-line tools are installed:
      - wget
      - tar
 3. Sufficient disk-space 
      - in your `/tmp` folder (>15GB) to build the container
      - in your working folder to store the container (~15GB)
      - for HippUnfold outputs (~1GB per subject) 
 4. Sufficient CPU and memory - the more you have, the faster it will run, but we recommend at least 4 CPU cores and 16GB memory.


## First time setup

Pull the container:

    apptainer pull jclauneurolab_autoafids_1.1.0.sif docker://dhananjhay/autoafids:1.1.0


Run AutoAFIDs without any arguments to print the short help:

    apptainer run -e jclauneurolab_autoafids_1.1.0.sif 

Use the `-h` option to get a detailed help listing:

    apptainer run -e jclauneurolab_autoafids_1.1.0.sif -h

Note that all the Snakemake command-line options are also available in
AutoAFIDs, and can be listed with `--help-snakemake`:

    apptainer run -e jclauneurolab_autoafids_1.1.0.sif --help-snakemake

Note: If you encounter any errors pulling the container from dockerhub, it may be because you are running 
out of disk space in your cache folders. Note, you can change these locations 
by setting environment variables, however, using a network file system for the folders may result in poor performance and/or errors e.g.:
    
    export APPTAINER_CACHEDIR=/YOURDIR/.cache/apptainer


## Running an example

Download and extract a single-subject BIDS dataset for this test:

    wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
    tar -xvf hippunfold_test_data.tar

This will create a `ds002168/` folder with a single subject, that has a 
both T1w and T2w images:

```
ds002168/
├── dataset_description.json
├── README.md
└── sub-1425
    └── anat
        ├── sub-1425_T1w.json
        ├── sub-1425_T1w.nii.gz
        ├── sub-1425_T2w.json
        └── sub-1425_T2w.nii.gz

2 directories, 6 files
```

Now let's run AutoAFIDs. 

    apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds002168 ds002168_autoafids participant -n --modality T1
Explanation:

Everything prior to the container (`jclauneurolab_autoafids_1.1.0.sif`) are arguments to apptainer, and after are to AutoAFIDs itself. The first three arguments to AutoAFIDs (as with any BIDS App) are the input
folder (`ds002168`), the output folder (`ds002168_autoafids`), and then the analysis level (`participant`). The `participant` analysis 
level is used in AutoAFIDs for performing any
participant-level processing. Here 
we used the T1w image. We also used the `--dry-run/-n`  option to 
just print out what would run, without actually running anything.


When you run the above command, a long listing will print out, describing all the rules that 
will be run. This is a long listing, and you can better appreciate it with the `less` tool. We can
also have the shell command used for each rule printed to screen using the `-p` Snakemake option:

    apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds002168 ds002168_autoafids participant -np | less


Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells AutoAFIDs how many cores to use.
 Using `--cores 8` means that AutoAFIDs will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (autoafids on a single subject) may take ~6 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.


    apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds002168 ds002168_autoafids participant -p --cores all


Note that you may need to adjust your [Singularity options](https://sylabs.io/guides/3.1/user-guide/cli/apptainer_run.html) to ensure the container can read and write to yout input and output directories, respectively. You can bind paths easily by setting an 
environment variable, e.g. if you have a `/project` folder that contains your data, you can add it to the `APPTAINER_BINDPATH` so it is available when you are running a container:

    export APPTAINER_BINDPATH=/data:/data



After this completes, you should have a `ds002168_autoafids` folder with outputs for the one subject.

## Exploring different options

If you alternatively want to run AutoAFIDs using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. The T2w image would be then first pre-processesed into T1w using SynthSR and then used as an input image for predicting AFIDs.

    apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds002168 ds002168_autoafids_t2w participant --modality T2w -p --cores all

Note that if you run with a different modality, you should use a separate output folder, since some of the files 
would be overwritten if not.



