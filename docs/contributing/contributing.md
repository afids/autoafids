# Contributing to AutoAFIDs

AutoAFIDs is a Python-based BIDS App for automatic anatomical fiducial detection. Development is managed using **pixi** for dependency and environment management.

> 🛠️ These instructions are intended for contributors making code changes to the AutoAFIDs codebase or using advanced features such as Snakemake cluster execution profiles.

---

## 📦 Prerequisites

Ensure **pixi** is installed on your system. You can install pixi by following the official guide:  
👉 [https://pixi.prefix.dev/latest/installation/](https://pixi.prefix.dev/latest/installation/)

Note: AutoAFIDs primarily supports T1-weighted images. For additional modalities (e.g., T2w), SynthSR will be triggered, though compatibility varies by operating system.

## 🧪 Setting Up the Development Environment

1. Clone the repository:

```bash
git clone https://github.com/afids/autoafids.git
cd autoafids
```

2. Create and activate the development environment:

```bash
pixi shell
```

3. Run AutoAFIDs with development dependencies:

```bash
./autoafids/run.py -h
```

---


## 🧹 Code Quality and Formatting

AutoAFIDs uses several tools to ensure clean and consistent code:

- [`isort`](https://github.com/pycqa/isort) – for sorting imports. 
- [`snakefmt`](https://github.com/snakemake/snakefmt) – for formatting Snakemake files
- [`black`](https://github.com/psf/black) – for formatting Python code

### Check formatting:

1. Activate the development environment:

```bash
pixi shell --environment dev
```

2. Run tools

```bash
isort autoafids/*.py -c
snakefmt autoafids --check
black autoafids --check
```

### Auto-fix formatting:

```bash
isort autoafids/*.py
snakefmt autoafids
black autoafids
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