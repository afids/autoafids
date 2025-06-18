# Running AutoAFIDs with Docker on Windows

> 🐳 These instructions assume Docker is already installed on your Windows system.  
> Docker also works on Linux and macOS with similar commands, but this guide assumes usage via the Windows Command Prompt (CMD).

---

## 🛠️ First-Time Setup

### 1. Open the Command Prompt and Pull Container

Open your Windows Command Prompt by clicking the bottom left `Windows` button and type `cmd` followed by `Enter`. This is where you will enter your AutoAFIDs commands. Feel free to make a new directory with `mkdir` or move to a directory you would like to work out of with `cd`, and for this example we will work from:

    cd c:\Users\dhananjhay\Downloads\

Pull the container (this will take some time and storage space, but like an installation it only needs to be done once and can then be run on many datasets):

    docker pull jclauneurolab/autoafids

Run AutoAFIDs without any arguments to print the short help:

    docker run -it --rm jclauneurolab/autoafids   

Use the `-h` option to get a detailed help listing:

    docker run -it --rm jclauneurolab/autoafids -h

Note that all the Snakemake command-line options are also available in
AutoAFIDs, and can be listed with `--help-snakemake`:

    docker run -it --rm jclauneurolab/autoafids --help-snakemake

### 2. Example

Download and extract a single-subject BIDS dataset for this test from [ds003653.tar](https://www.dropbox.com/scl/fi/phmmofiy4q6o1k01rs6c4/ds003653.tar?rlkey=bpa8fxfl0lyrdc38fs6aowta7&st=zvhpqsga&dl=1"). Here we will also assume you chose to save and extract to the directory `c:\Users\dhananjhay\Downloads\`.

This contains a `ds003653/` directory with a single subject, that has a both T1w and T2w images. 

```
ds003653/
├── dataset_description.json
├── README
└── sub-718211
    └── ses-01
        ├── anat
        │   ├── sub-718211_ses-01_T1w.json
        │   ├── sub-718211_ses-01_T1w.nii.gz
        │   ├── sub-718211_ses-01_T2w.json
        │   └── sub-718211_ses-01_T2w.nii.gz
        ├── sub-718211_ses-01_scans.json
        └── sub-718211_ses-01_scans.tsv

3 directories, 8 files
```

Now let's run AutoAFIDs on the test dataset. Docker will need read/write access to the input and output directories, respectively. This is achieved with the `-v` flag. This 'binds' or 'mounts' a directory to a new directory inside the container.

    docker run -it --rm -v c:\Users\dhananjhay\Downloads\ds003653:/bids -v c:\Users\dhananjhay\Downloads\ds003653_autoafids:/output jclauneurolab/autoafids /bids /output participant -n

Explanation: 

`-v c:\Users\dhananjhay\Downloads\ds003653:/bids` tells Docker to mount the directory `c:\Users\dhananjhay\Downloads\ds003653` into a new directory inside the container named `/bids`. We then do the same things for our output directory named `ds003653_autoafids`, which we mount to `/output` inside the container. These arguments are not specific to AutoAFIDs but rather are general ways to use Docker. You may want to familiarize yourself with [Docker options](https://docs.docker.com/engine/reference/run/).

Everything after we specified the container (`jclauneurolab/autoafids`) are arguments to AutoAFIDs itself. The first of these arguments (as with any BIDS App) are the input directory (`/bids`), the output directory (`/output`), and then the analysis level (`participant`). The `participant` analysis 
level is used in AutoAFIDs for performing any participant-level processing. We also used the `--dry-run/-n`  option to just print out what would run, without actually running anything.

When you run the above command, a long listing will print out, describing all the rules that 
will be run. Now, to actually run the workflow, we need to specify how many cores to use and leave out
the dry-run option.  The Snakemake `--cores` option tells AutoAFIDs how many cores to use.
 Using `--cores 8` means that AutoAFIDs will only make use of 8 cores at most. Generally speaking 
you should use `--cores all`,  so it can make maximal use of all the CPU cores it has access to on your system. This is especially 
useful if you are running multiple subjects. 

Running the following command (autoafids on a single subject) may take ~6 minutes if you have 8 cores, shorter if you have more 
cores, but could be much longer (several hours) if you only have a single core.

    docker run -it --rm -v c:\Users\dhananjhay\Downloads\ds003653:/bids -v c:\Users\dhananjhay\Downloads\ds003653_autoafids:/output jclauneurolab/autoafids /bids /output participant -p --cores all

After this completes, you should have a `ds003653_autoafids` directory with outputs for the one subject.

## Exploring different options

If you alternatively want to run AutoAFIDs using a different modality, e.g. the high-resolution T2w image
in the BIDS test dataset, you can use the `--modality T2w` option. The T2w image would be then first pre-processesed into T1w using SynthSR and then used as an input image for predicting AFIDs.

    docker run -it --rm -v c:\Users\dhananjhay\Downloads\ds003653:/bids -v c:\Users\dhananjhay\Downloads\ds003653_autoafids_t2w:/output jclauneurolab/autoafids /bids /output participant --modality T2w -p --cores all

Note that if you run with a different modality, you should use a separate output directory, since some of the files 
would be overwritten if not.


