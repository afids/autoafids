# Changelog

See [semantic versioning](https://semver.org/spec/v2.0.0.html) for the rationale behind
the version numbers.


## [1.1.0] (2025-03-07)

### Added

- Docker support for easier deployment.
- Restricted supported Python versions to improve compatibility.
- Upgraded Snakebids to the latest version.
- Reorganized outputs by bifurcating work and root for better navigation.
- Introduced stereotaxy support for precise coordinate-based processing.
- Enhanced multi-modality support by adding SynthSR for non-T1w modalities (T2w, FLAIR).
- Renamed --profile to --procprofile to prevent conflicts when running AutoAFIDS on compute clusters.

### Changed

- Removed model extraction on a per-subject basis to optimize performance.

### Fixed

- Resolved incorrect use of the expand function for better processing.
- Eliminated hardcoded T1w modality to support flexible pipeline configurations.
- Fixed minor pipeline issues to improve stability and performance.


## [1.0.0] (2024-12-06)

### Added

- Integrated preprocessing pipeline with apply workflow
- Uploaded models on cloud storage
- Optimized model loading and configuration
