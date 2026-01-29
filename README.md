SIFS: Spatially Informed Feature Selection for MALDI-MSI

![SIFS](SIFS_logo.png)


**SIFS** is an R package implementing *Spatially Informed Feature Selection* for mass spectrometry imaging (MSI), designed to prioritize **m/z features** that are informative **in their spatial context** (e.g., co-localization with neuropathology annotations and spatial molecular patterns), rather than relying on intensity differences alone.

This repository accompanies the manuscript from the CeMOS Mannheim group and provides a reproducible implementation of the SIFS workflow used for spatially aware **m/z feature selection** (e.g., reducing high-dimensional MSI spectra to a compact, informative set),


---

## Key ideas

SIFS is built around the following concepts:

- **Spatial context matters:** MSI signals are evaluated not only by magnitude but by their *spatial organization*.
- **Co-localization with annotations:** Features can be prioritized based on agreement with histology/neuropathology annotations/labels (where provided).
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

* matrix/statistics utilities (packages: MALDIquant (1.22.3), MALDIquantForeign (0.14.1)),
* spatial/image operations (packages: moleculaR (0.9.5), spatstat (3.5-0), spatstat.geom (3.1-9)),
* machine learning wrappers (optional; SIFS itself focuses on feature selection) 
 * package requirements for ML inside Rstudio: (reticulate (1.44.1), python (3.12), SHAP (0.5) xgboost (3.1.3)).
* 

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

# 1) Run spatially informed feature selection: 
# selected <- HMCS_calculatorBinary(dsc_df, focuROI = "VT")
## Print(head(selected, 5))
## Note: returns a dataframe of three columns: mzList, dsc_mpm, HMCS_VT
## this dataframe must be ordered based on the last column `(HMCS_VT)` then it can be fed to the downstream ML/DL or benchmarking pipelines.  

# 2) Subset the feature matrix
# X_red <- X[, selected$mz_index]

# 3) Use the reduced features in your classifier of choice
# model <- train_classifier(X_red, y)  # e.g., xgboost / random forest / etc.
```
---

## Citation

If you use SIFS in academic work, please cite the associated manuscript:

> Citation will be included once the article is published. 


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
