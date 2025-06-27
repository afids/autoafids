# Running AutoAFIDs with Apptainer (Singularity)

## Prerequisites

Before running AutoAFIDs with Apptainer (formerly Singularity), ensure the following requirements are met:

1. **Apptainer is installed.**  
   See the official [Apptainer installation guide](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages) for instructions.

2. **Required command-line tools are installed:**
   - `wget`
   - `tar`

3. **Sufficient disk space is available:**
   - At least **15 GB** free in `/tmp` to build the container
   - At least **15 GB** in your working directory to store the `.sif` container
   - Approximately **1 GB per subject** for AutoAFIDs output

4. **Adequate CPU and memory:**
   - Minimum: 4 CPU cores and 16 GB RAM
   - More resources will improve performance


## First-Time Setup

### 🔄 Pull the container

```bash
apptainer pull autoafids.sif docker://jclauneurolab/autoafids
```

### 🚀 Run AutoAFIDs

Run without arguments to see a short help message:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif
```

Get detailed help with the `-h` flag:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif -h
```

List available Snakemake-specific options using:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif --help-snakemake
```

---

### ⚠️ Troubleshooting: Pull Errors

If you encounter errors while pulling the container, it may be due to insufficient disk space in your Apptainer cache directories.

You can redirect the cache location by setting the following environment variable:

```bash
export APPTAINER_CACHEDIR=/your/custom/path/.cache/apptainer
```

> 💡 **Tip:** Avoid using network file systems (e.g., NFS) for the cache directory, as this may lead to performance issues or errors during image creation.



## Running an Example

### 📥 Download a Sample Dataset

Download and extract a single-subject BIDS dataset:

```bash
wget "https://www.dropbox.com/scl/fi/phmmofiy4q6o1k01rs6c4/ds003653.tar?rlkey=bpa8fxfl0lyrdc38fs6aowta7&st=zvhpqsga&dl=1" -O ds003653.tar
tar -xvf ds003653.tar
```

This will create a `ds003653/` folder containing a single subject with both T1w and T2w images:

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
```

---

### 🚀 Run AutoAFIDs

Run the following command to perform a dry-run with the T1w image:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds003653 ds003653_autoafids participant -n --modality T1
```

#### 🔍 Explanation

Everything **before** the container name is passed to Apptainer. Everything **after** is passed to AutoAFIDs:

- `ds003653`: input folder
- `ds003653_autoafids`: output folder
- `participant`: analysis level
- `--modality T1`: use T1w image
- `-n`: dry-run mode (no processing, just show steps)

To inspect the rule graph more interactively:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds003653 ds003653_autoafids participant -np | less
```

To actually run the pipeline, remove `-n` and specify the number of cores:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds003653 ds003653_autoafids participant -p --cores all
```

> 💡 We recommend `--cores all` to utilize all available CPU cores. This speeds up processing, especially for multi-subject datasets.

---

### ⚙️ Apptainer Bind Paths

Make sure input/output folders are accessible inside the container. You can bind host paths using:

```bash
export APPTAINER_BINDPATH=/data:/data
```

Adjust paths as needed for your system.

---

### 📂 After Completion

Upon completion, the `ds003653_autoafids/` folder will contain output files for the processed subject.

---

## Exploring Different Options

To run AutoAFIDs using the T2w image (or CT and FLAIR) instead:

```bash
apptainer run -e jclauneurolab_autoafids_1.1.0.sif ds003653 ds003653_autoafids_t2w participant --modality T2w -p --cores all
```

> ⚠️ Use a **separate output folder** when running with a different modality, as outputs may be overwritten.

When using `--modality T2w`, the T2w image will first be synthesized into a T1-like image using **SynthSR**, then used for fiducial prediction.




