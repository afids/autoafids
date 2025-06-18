# Installation

## AutoAFIDs: A BIDS App for Automatic Landmark Localization

***AutoAFIDs*** is a BIDS App that automatically detects anatomical landmarks in the human brain. It is trained to predict the coordinates of a standardized set of well-validated anatomical landmarks (e.g., the anterior and posterior commissures).

---

## 🧰 System Requirements

we support three modes for installing ***AutoAFIDs***:

| Method                | OS Support                  | Container Needed | Editable | Cluster Profiles | Notes                                |
|-----------------------|-----------------------------|------------------|----------|------------------|---------------------------------------|
| **Conda**             | Linux, macOS (Intel & M2 only)   | ❌ No            | ✅ Yes   | ✅ Yes           | Recommended for most local setups     |
| **Docker**            | Windows, macOS (Intel/Silicon) | ✅ Yes     | ❌ No    | ❌ No            | Works cross-platform (slow on M1)     |
| **Singularity/Apptainer** | Linux only             | ✅ Yes           | ❌ No    | ❌ No            | Good for HPC or shared environments   |

> 💡 **Note:** GPU is *not* required.  
> ⚠️ **Apple Silicon (M1):** Native Conda support is not available. You can use Docker with emulated `amd64`, but performance will be slow.

If you're interested in contributing to development, see the [Contributing Guide](https://autoafids.readthedocs.io/en/latest/contributing/contributing.html).

---

## 📂 Inputs and Outputs

- **Input:** BIDS-formatted dataset containing structural MRI images (best on T1-weighted modalities).
- **Output:** Automatically placed anatomical fiducials (AFIDs) saved as structured coordinate files and visualizations (if using `--fidqc` flag).

---

## 🚀 Installation Methods

### 🔧 1. Conda Environment (Linux/macOS)

AutoAFIDs can be installed via [Conda](https://docs.conda.io/) for container-free execution on Linux and Intel macOS systems.

📘 [Conda Installtion Instructions](https://autoafids.readthedocs.io/en/latest/getting_started/conda.html)


**✅ Pros:**
- Easy to install and run
- Fully editable code
- Works with Snakemake profiles

**⚠️ Cons:**
- Not supported on Windows or Apple Silicon (M1) system
- Requires Conda or Miniconda

---

### 🐳 2. Docker Container (Windows/Linux/macOS Intel)

A pre-built Docker container is available with all dependencies bundled.

📘 [Docker Setup Instructions](https://autoafids.readthedocs.io/en/latest/getting_started/docker.html)

**✅ Pros:**
- Cross-platform support, including limited M1 compatibility via emulation
- All dependencies packaged in a single file

**⚠️ Cons:**
- Cannot edit code inside the container
- Not compatible with Snakemake cluster execution
- May not be usable on shared systems (e.g., university clusters)

---

### 🧪 3. Singularity / Apptainer (Linux)

Use [Singularity](https://apptainer.org/) (now Apptainer) to convert the Docker image into a `.sif` file for HPC or institutional clusters.

📘 [Singularity Setup Instructions](https://autoafids.readthedocs.io/en/latest/getting_started/singularity.html)

**✅ Pros:**
- Ideal for shared Linux systems
- All dependencies packaged in a single file

**⚠️ Cons:**
- Cannot edit code inside the container
- Not compatible with Snakemake cluster execution
