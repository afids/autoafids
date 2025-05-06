# Running AutoAFIDS with Conda

AutoAFIDs can be installed and run using Conda on **Linux and macOS** systems.

**Note:** Conda installation is **not supported on Windows** at this time. If you are on Windows, please refer to the [Docker instructions](docker.md) instead.

---

## For Users: Installing AutoAFIDs via Conda

These steps are intended for **end users** who simply want to run AutoAFIDs and get their AFIDs.

### 1. Install Conda (if not already installed)

Follow the instructions at the official Conda documentation:
[https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

---

### 2. Create and activate a new Conda environment

```bash
conda install mamba -c conda-forge
mamba create --name autoafids-env -c khanlab -c conda-forge -c bioconda autoafids
mamba activate autoafids-env
```

---

### 3. Test the installation

Run the following command to verify the installation:

```bash
autoafids -h
```

You should see a help message listing all available command-line options.  
If this runs successfully, you’re ready to start processing data with AutoAFIDs!

## Running an example

You can try AutoAFIDs on a sample dataset to make sure everything works as expected.

First, download and extract a single-subject BIDS dataset for this test:

```bash
wget https://www.dropbox.com/s/mdbmpmmq6fi8sk0/hippunfold_test_data.tar 
tar -xvf hippunfold_test_data.tar
```

This will create a `ds002168/` folder with a single subject that has both T1w and T2w images:

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

### Run the full AutoAFIDs BIDS pipeline

Running AutoAFIDs:

```bash
autoafids ds002168 ds002168_autoafids participant --cores all --use-conda
```

This should run the full pipeline and place results in a new `ds002168_autoafids/` folder.


## Cache Directory

When running, AutoAFIDs automatically downloads and caches necessary SynthSR repo and CNN model to speed up subsequent runs.

By default, these are stored in the following directory:

```bash
~/.cache/autoafids/
```

You can override this default cache location by setting the `AUTOAFIDS_CACHE_DIR` environment variable:

```bash
export AUTOAFIDS_CACHE_DIR=/path/to/custom/cache
```

This is useful when working on shared systems, when home directory storage is limited, or if you wish to isolate data per project or user.

## For Developers & Contributors

These steps are intended for people who want to contribute to the development of AutoAFIDs or explore its internals.

### 1. Clone the AutoAFIDs GitHub repository

```bash
git clone https://github.com/afids/autoafids.git
cd autoafids
```

---

### 2. Create and activate a new Conda environment

```bash
conda env create -f autoafids/autoafids-dev.yml
conda activate autoafids-dev
```

---

### 3. Run the development version of AutoAFIDs

You can run AutoAFIDs directly from the source directory using:

```bash
./autoafids/run.py -h
```

This should print out the available command-line arguments for the tool.  
You’re now set up for development and contribution!

---

## Troubleshooting

If you encounter issues while setting up AutoAFIDs via Conda:

- Make sure you’re using the latest Conda:
  ```bash
  conda update -n base -c defaults conda
  ```
- Double-check that your environment is activated (`conda activate autoafids-env` or `autoafids-dev`)
- Try creating a fresh environment if problems persist
- Search for similar issues or open a new one in the [GitHub issues](https://github.com/afids/autoafids/issues) page

---

Happy AFIDs placing!
