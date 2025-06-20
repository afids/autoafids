# Automatic Anatomical Fiducials (AutoAFIDs)

[![Documentation Status](https://readthedocs.org/projects/autoafids/badge/?version=latest)](https://autoafids.readthedocs.io/en/stable/?badge=stable)
![Version](https://img.shields.io/github/v/tag/afids/autoafids?label=version)
![Python3](https://img.shields.io/badge/python-_3.9_|_3.10_|_3.11_|_3.12-blue.svg)
[![Tests](https://github.com/afids/autoafids/actions/workflows/lint-and-dryrun-testing.yml/badge.svg?branch=main)](https://github.com/afids/autoafids/actions/workflows/lint-and-dryrun-testing.yml?query=branch%3Amain)
![Docker Pulls](https://img.shields.io/docker/pulls/jclauneurolab/autoafids)

Developed by the AIMS Lab at the Robarts Research Institute  
*2023–2025*

> ⚠️ This package is under active development. While stable and reproducible, users are encouraged to report any bugs or unexpected behavior to the development team.

---

## 📑 Table of Contents

- [🧠 What is AutoAFIDs?](#what-is-autoafids)
- [⚙️ Workflow Overview](#workflow-overview)
- [🔍 Known Issues](#known-issues)
- [📖 Full Documentation](#full-documentation)
- [💬 Questions, Issues, and Feedback](#questions-issues-and-feedback)
- [📚 Relevant Papers](#relevant-papers)

---
## What is AutoAFIDs?

**AutoAFIDs** is a BIDS App for automatic detection of anatomical fiducials (AFIDs) on MRI scans. We make use of these AFIDs in various neuroimaging applications such as image registration, quality control, and neurosurgical targeting.

AutoAFIDs leverages:

- A 3D [U-Net architecture](https://arxiv.org/abs/1505.04597) for landmark localization
- The [Snakemake](https://snakemake.readthedocs.io/) and [SnakeBIDS](https://snakebids.readthedocs.io/en/stable/) workflow frameworks for reproducible processing
- The AFIDs protocol developed by the [AFIDs community](https://github.com/afids)

The software is modality-aware, but best supports T1-weighted MRI scans. It also includes QC visualizations for quality assurance.

---

## Workflow Overview

Below is a high-level summary of the AutoAFIDs processing pipeline:

![Pipeline Overview](https://raw.githubusercontent.com//afids/autoafids/master/docs/images/dag.png)

1. Preprocess BIDS input files based on image modality (T1w, T2w)
2. Load trained fiducial models (32 AFID points)
3. Run patch-based U-Net inference per AFID
4. Generate predictions

---

## Known Issues

- `gen_fcsv` rule is currently sequential and may benefit from AFID-level parallelization
- T1w-like scan synthesis is available (via SynthSR [Iglesias et al., 2023](https://www-science-org.proxy1.lib.uwo.ca/doi/10.1126/sciadv.add3607)) but requires millimetric validation and may not work for all operating systems

---

## Full Documentation

👉 [autoafids.readthedocs.io](https://autoafids.readthedocs.io/en/)

Includes installation instructions, usage examples, and advanced configuration.

---

## Questions, Issues, and Feedback

- Open an issue on the [GitHub issues page](https://github.com/afids/autoafids/issues)
- Start a conversation on the [Discussions board](https://github.com/afids/autoafids/discussions)
- Contact: [ataha24@uwo.ca](mailto:ataha24@uwo.ca) or [dbansal7@uwo.ca](mailto:dbansal7@uwo.ca)

We welcome feedback, contributions, and collaboration!


## Relevant Papers

AutoAFIDs builds upon a series of foundational works that introduce, validate, and apply the anatomical fiducials (AFIDs) framework for neuroimaging quality control, registration evaluation, and surgical planning.

---

### Methodological Foundations

- **Lau et al., 2019**  
  *A framework for evaluating correspondence between brain images using anatomical fiducials.*  
  *Human Brain Mapping*, 40(14), 4163–4179.  
  [DOI: 10.1002/hbm.24695](https://doi.org/10.1002/hbm.24695)

---

### Clinical Applications

- **Abbass et al., 2022**  
  *Application of the anatomical fiducials framework to a clinical dataset of patients with Parkinson’s disease.*  
  *Brain Structure & Function*, 227(1), 393–405.  
  [DOI: 10.1007/s00429-021-02408-3](https://doi-org.proxy1.lib.uwo.ca/10.1007/s00429-021-02408-3)

- **Taha et al., 2022**  
  *An indirect deep brain stimulation targeting tool using salient anatomical fiducials.*  
  *Neuromodulation: Journal of the International Neuromodulation Society*, 25(8), S6–S7.

---

### Dataset & Resource Publication

- **Taha et al., 2023**  
  *Magnetic resonance imaging datasets with anatomical fiducials for quality control and registration.*  
  *Scientific Data*, 10(1), 449.  
  [DOI: 10.1038/s41597-023-02330-9](https://doi.org/10.1038/s41597-023-02330-9)

---

### Registration and Localization Accuracy

- **Abbass et al., 2025**  
  *The impact of localization and registration accuracy on estimates of deep brain stimulation electrode position in stereotactic space.*  
  *Imaging Neuroscience.*  
  [DOI: 10.1162/imag_a_00579](https://doi.org/10.1162/imag_a_00579)

---

> 📌 If you use AutoAFIDs or AFIDs data in your research, please cite the relevant papers above.