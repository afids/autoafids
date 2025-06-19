# Contributing to AutoAFIDs

AutoAFIDs is a Python-based BIDS App for automatic anatomical fiducial detection. Development is managed using **Conda** for dependency and environment management.

> 🛠️ These instructions are intended for contributors making code changes to the AutoAFIDs codebase or using advanced features such as Snakemake cluster execution profiles.

---

## 📦 Prerequisites

Ensure **Conda** is installed on your system. You can install Miniconda or Anaconda by following the official guide:  
👉 [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

<<<<<<< Updated upstream
Note: AutoAFIDs primarily supports T1-weighted images. For additional modalities (e.g., T2w), `SynthSR` will be triggered, though compatibility varies by operating system.

=======
>>>>>>> Stashed changes
---

## 🧪 Setting Up the Development Environment

1. Clone the repository:

```bash
git clone https://github.com/afids/autoafids.git
cd autoafids
```

2. Create and activate the development environment:

```bash
mamba env create -f autoafids-dev.yml
eval "$(mamba shell hook --shell zsh)"
mamba activate autoafids-dev
```

3. Run AutoAFIDs with development dependencies:

```bash
./autoafids/run.py -h
```

---


## 🧹 Code Quality and Formatting

AutoAFIDs uses several tools to ensure clean and consistent code:

- [`ruff`](https://github.com/charliermarsh/ruff) – for linting and formatting
- [`snakefmt`](https://github.com/snakemake/snakefmt) – for formatting Snakemake files
- [`yamlfix`](https://github.com/lyz-code/yamlfix) – for YAML cleanup

### Check formatting:

```bash
ruff check .
snakefmt --check Snakefile
yamlfix --check config/
```

### Auto-fix formatting:

```bash
ruff check . --fix
snakefmt Snakefile
yamlfix config/
```

---

## 🧪 Dry Run / Workflow Testing

To test the Snakemake workflow without running any jobs, use the built-in **dry-run** feature:

```bash
./autoafids/run.py tests/data tests/output participant -n
```

The `tests/data` directory contains a lightweight fake BIDS dataset (with zero-byte files) useful for testing the pipeline logic.

You can simulate more complex scenarios by modifying the `tests/` configuration and rerunning dry-runs with `--modality` and other options.

---

## 🙋 Questions, Issues, and Feedback

We welcome all forms of feedback and contributions.

- 💬 For questions or suggestions, use the [Discussions](https://github.com/afids/autoafids/discussions) page.
- 🐛 For bugs or feature requests, open an issue on the [Issues](https://github.com/afids/autoafids/issues) page.
- 📧 You may also contact [dbansal7@uwo.ca](mailto:dbansal7@uwo.ca) or reach out to **Alaa Taha**, the Science Lead.

Thanks for helping improve AutoAFIDs!