# Running AutoAFIDs with pixi

AutoAFIDs can be installed and run using pixi on **Linux** system only. Pixi will manage all Python dependencies and non-python dependencies (c3d, greedy, ANTS) through conda environments.

**Note:** Pixi installation is **not supported on Windows** at this time. If you are on Windows, please refer to the [Docker instructions](docker.md) instead.


## For Users: Installing AutoAFIDs via Pixi


### Steps
 1. Install pixi (if not already installed):
    ```bash
    curl -fsSL https://pixi.sh/install.sh | bash
    ```

 2. Clone the repository and install dependencies:
    ```bash
    git clone https://github.com/afids/autoafids.git
    cd autoafids
    pixi install
    ```

## Development
For development work, use the development environment which includes additional tools like formatters and linters:

```bash
pixi install --environment dev
```

Quality checks can be run using:
```bash
pixi run --environment dev quality_check  # Check code formatting
pixi run  --environment dev quality_fix    # Fix code formatting
```   

## Usage

### Test the installation

Run the following command to verify the installation:

```bash
pixi run autoafids -h
```

You should see a help message listing all available command-line options.  
If this runs successfully, you’re ready to start processing data with AutoAFIDs!

## Running an example

You can try AutoAFIDs on a sample dataset to make sure everything works as expected.

First, download and extract a single-subject BIDS dataset for this test:

```bash
wget "https://www.dropbox.com/scl/fi/phmmofiy4q6o1k01rs6c4/ds003653.tar?rlkey=bpa8fxfl0lyrdc38fs6aowta7&st=zvhpqsga&dl=1" -O ds003653.tar
tar -xvf ds003653.tar
```

This will create a `ds003653/` folder with a single subject that has both T1w and T2w images:

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

### Run the full AutoAFIDs BIDS pipeline

By default (Linux or Intel-based macOS), run:

```bash
pixi run autoafids ds003653 ds003653_autoafids participant --cores all
```

This will run the full pipeline and place results in a new `ds003653_autoafids/` folder.


## Cache Directory

When running, AutoAFIDs automatically downloads and caches the necessary CNN model to speed up subsequent runs.

By default,  it's stored in the following directory:

```bash
~/.cache/autoafids/
```

You can override this default cache location by setting the `AUTOAFIDS_CACHE_DIR` environment variable:

```bash
export AUTOAFIDS_CACHE_DIR=/path/to/custom/cache
```

This is useful when working on shared systems, when home directory storage is limited, or if you wish to isolate data per project or user.

Happy AFIDs placing!