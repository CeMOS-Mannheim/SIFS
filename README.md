# SIFS: Spatially Informed Feature Selection for MALDI-MSI

**SIFS** is an R package implementing *Spatially Informed Feature Selection* for mass spectrometry imaging (MSI), designed to prioritize **m/z features** that are informative **in their spatial context** (e.g., co-localization with neuropathology annotations and spatial molecular patterns), rather than relying on intensity differences alone.

This repository accompanies the manuscript from the CeMOS Mannheim group and provides a reproducible implementation of the SIFS workflow used for spatially aware **m/z feature selection** (e.g., reducing high-dimensional MSI spectra to a compact, informative set),

---

## Key ideas

SIFS is built around the following concepts:

- **Spatial context matters:** MSI signals are evaluated not only by magnitude but by their *spatial organization*.
- **Co-localization with annotations:** Features can be prioritized based on agreement with histology/neuropathology labels (where provided).
- **Compact, performant feature sets:** SIFS supports aggressive feature reduction (e.g., to 256 m/z) for faster learning and improved generalization in many models.

---

## Installation

### Stable installation from GitHub

Install the development version directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("CeMOS-Mannheim/SIFS")
````

Load the package:

```r
library(SIFS)
```

---

## Dependencies

SIFS is an R package and depends on common scientific R libraries. The exact dependency list is recorded in `DESCRIPTION`, but you can typically expect packages for:

* matrix/statistics utilities,
* spatial/image operations,
* machine learning wrappers (optional; SIFS itself focuses on feature selection).

If you encounter installation issues on Linux/macOS related to system libraries, please ensure you have a working C/C++ toolchain and standard build tools for R packages.

---

## Quick start (minimal example)

Below is a schematic example showing the typical steps. The exact function names may differ depending on the finalized exported API; this README will be updated once the public interface is frozen.

```r
library(SIFS)

# Example inputs (conceptual)
# X: matrix [n_pixels x n_mz] intensity matrix
# coords: matrix/data.frame [n_pixels x 2] of spatial coordinates
# y: optional vector [n_pixels] with neuropathology labels / ROI membership

# 1) Run spatially informed feature selection
# selected <- sifs_select(X, coords = coords, y = y, k = 256)

# 2) Subset the feature matrix
# X_red <- X[, selected$mz_index]

# 3) Use the reduced features in your classifier of choice
# model <- train_classifier(X_red, y)  # e.g., xgboost / random forest / etc.
```
---

## Citation

If you use SIFS in academic work, please cite the associated manuscript:

> Citation will be included once paper is published. 


---

## License

This project is distributed under the license specified in `LICENSE` (and `DESCRIPTION`).

---

## Contributing

Issues and pull requests are welcome.

If you report a bug, please include:

* a minimal reproducible example,
* `sessionInfo()`,
* a short description of the MSI data structure (dimensions, file format, preprocessing steps).

---

## Contact

CeMOS Mannheim / Mannheim research team (see repository maintainers).

---
