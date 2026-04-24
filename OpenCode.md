# LLM Wiki — OpenCode · SIFS

> A personal knowledge base maintained by **OpenCode**.  
> Based on Andrej Karpathy's LLM Wiki pattern.  
> Domain: **SIFS** — Spatially Informed Feature Selection for Mass Spectrometry Imaging (MSI).  
> Author: Shad Arif Mohammed · CeMOS, Mannheim University of Applied Sciences.

---

## 1 · Purpose

This wiki is a structured, interlinked knowledge base for the **SIFS** R package project. It supports active development of vignettes, package functions, and documentation — enabling efficient LLM-assisted editing with the `graphify` library to minimise token usage.

**OpenCode** maintains the wiki.  
**The human** curates the source code, asks questions, directs edits, and decides what ships to GitHub.

> ⚠ **GitHub boundary:** All wiki files (`wiki/`), raw working data (`data/`, `raw/`), interim outputs (`output/`), and `graphify-out/` are listed in `.gitignore` and must **never** be pushed to the public repository.

### Primary outputs

| Output | Format | Location |
|--------|--------|----------|
| Vignette Rmd notebooks | `.Rmd` | `vignettes/` |
| Exported package data | `.rda` | `data/` |
| Raw intermediate objects | `.rds`, `.csv` | `raw/` / `data/output/` |
| Wiki concept pages | Markdown | `wiki/` |
| Graphify dependency graphs | JSON / PNG | `graphify-out/` |

---

## 2 · Project Directory Structure

```
SIFS/
├── R/                    ← exported package functions (edit here, ships to GitHub)
├── man/                  ← Roxygen-generated documentation (auto-generated, ships)
├── vignettes/            ← Rmd walkthroughs (ships to GitHub)
│   └── 01_SIFS_workflow1.Rmd
├── data/                 ← LazyData .rda objects (gitignored)
│   ├── input/
│   │   ├── 0_mis_files/
│   │   ├── 01_spmatObjects/
│   │   │   ├── trainset_tissue/
│   │   │   └── testset_tissue/
│   │   ├── 04_y_lbls/
│   │   └── 06_fwhm_object/
│   └── output/
│       ├── 01_pos/
│       ├── 02_neg/
│       ├── models/
│       ├── SHAP_outputs/
│       │   ├── VT_SHAPs/
│       │   ├── Nec_SHAPs/
│       │   └── PreNec_SHAPs/
│       ├── SPPs_MPMs/
│       └── plots/
├── raw/                  ← immutable source objects (gitignored)
├── wiki/                 ← OpenCode knowledge pages (gitignored)
│   ├── index.md
│   └── log.md
├── graphify-out/         ← graphify artefacts (gitignored)
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── SIFS.Rproj
└── .gitignore
```

> **Rule:** Never modify files inside `raw/`. All generated outputs go to `data/output/`. Only `R/`, `man/`, `vignettes/`, `DESCRIPTION`, `NAMESPACE`, and `README.md` ship to GitHub.

---

## 3 · Package Overview

### What SIFS does

SIFS implements a two-stage spatially informed feature selection pipeline for MSI data, prioritising m/z features that are meaningful in their *spatial* context rather than purely by statistical discriminability.

**Stage 1 — SHAP-based feature ranking (Python/XGBoost):**
- An XGBoost classifier is trained on a sparse intensity matrix (`spmat`) from a training-set tissue.
- SHAP (SHapley Additive exPlanations) values are extracted on a held-out testing-set tissue.
- Features with zero SHAP values across all ROI classes are removed via `Zero_remover()`.
- Remaining features are assigned to an ROI class (e.g. VT, Nec, PreNec) by their maximum SHAP value.

**Stage 2 — HMCS spatial scoring (R/moleculaR):**
- For each SHAP-selected m/z, a Spatial Peak Profile (SPP) is computed via `moleculaR::searchAnalyte()`.
- A Molecular Probability Map (MPM) / hotspot mask is derived via `moleculaR::probMap()`.
- The **Hotspot-Map Co-localisation Score (HMCS)** is calculated as the Dice Similarity Coefficient (DSC) between each m/z hotspot map and each annotated tissue ROI window.
- `HMCS_calculatorBinary()` scores features by `DSC(mpm, VT) − DSC(mpm, Nec)` to yield a ranked list per ROI.
- The top HMCS-ranked features are the final SIFS output — spatially validated, transferable m/z features.

### Key concepts

| Term | Meaning |
|------|---------|
| `spmat` | Sparse intensity matrix: named list with `$spmat` (dgCMatrix), `$mzAxis`, `$coordinates` |
| `mis` | Measurement Information Structure: spatial metadata + tissue boundary + ROI annotations |
| `spwin` | Spatial window (`moleculaR::createSpatialWindow`): bounding window of pixel coordinates |
| `SPP` | Spatial Peak Profile: ion image for a specific m/z with Gaussian weighting |
| `MPM` | Molecular Probability Map: kernel-density hotspot probability surface |
| `DSC` | Dice Similarity Coefficient: spatial overlap score ∈ [0, 1] |
| `HMCS` | Hotspot-Map Co-localisation Score: signed DSC contrast between ROIs |
| `FWHM` | Full-Width at Half Maximum: instrument peak width used for m/z search tolerance |
| `SHAP` | SHapley Additive exPlanations: feature importance from XGBoost |
| ROI | Region of Interest — VT (Viable Tumour), Nec (Necrosis), PreNec (Pre-Necrotic) |

### Tissue sample reference

| Ref ID | Tissue label | Role |
|--------|-------------|------|
| `t3` | XIII-t3 83318 | **Training set** |
| `t4` | XV-t1 82558 | **Testing set** |

---

## 4 · Exported R Functions

Functions live in `R/`. Document with Roxygen2; regenerate `man/` via `devtools::document()`.

| Function | Description |
|----------|-------------|
| `Zero_remover(df, col_idx)` | Removes rows where all selected columns are zero; used for SHAP filtering |
| `HMCS_calculatorBinary(DSC_df, focusROI)` | Computes binary HMCS score: DSC(focus ROI) − DSC(contrast ROI) |

> When adding new exported functions, update this table and create a wiki page under `wiki/fn-<name>.md`.

---

## 5 · Vignette Structure — `01_SIFS_workflow1.Rmd`

The primary vignette walks through a complete end-to-end SIFS run on two GBM tissue samples. Sections in order:

1. **Install & load** — `remotes::install_github("CeMOS-Mannheim/SIFS")`
2. **Global variables** — `root_dir`, `Alpha`, `current_date`, `ionisation`, `ROIs`, `ROI`, `num_comparison`
3. **I/O paths** — `input_main`, `output_main` (switched on ionisation mode `pos`/`neg`)
4. **SHAP mode selection** — `sel_mzModes`, `SelSHAP` (PosSHAP / NegSHAP)
5. **Tissue cohort metadata** — `t_names`, `trainSetNames`, `testSetNames`
6. **Load MIS files** — `t3_mis`, `t4_mis` from `.rds`
7. **Set ROI names** — `TrainingSetROI_name = "83318"`, `TestingSetROI_name = "82558"`
8. **Load SPMAT objects** — `t3_spmat_pos`, `t4_spmat_pos`
9. **Load labels** — `lbls_pos_ALL`, split by tissue, ROI-conditional assignment of `t3_lbls` / `t4_lbls`
10. **Python / SHAP block** — XGBoost training with optional hyperparameter tuning, SHAP extraction, save CSVs to `output/SHAP_outputs/`
11. **SHAP filter** — `Zero_remover()`, feature class assignment by max SHAP
12. **Spatial windows & FWHM** — `moleculaR::createSpatialWindow()`, load `fwhm_pos`
13. **Select mzList** — choose m/z list for the current ROI from SHAP-assigned lists
14. **SPP + MPM loop** — `moleculaR::searchAnalyte()` + `moleculaR::probMap()` for `t3` and `t4`
15. **Save SPP/MPM** — `saveRDS()` to `output/SPPs_MPMs/`
16. **DSC calculation** — `moleculaR::dsc()` for each MPM against each ROI window
17. **HMCS scoring** — `HMCS_calculatorBinary()`, produce ranked `dsc_df`
18. **Visualisation** — `knitr::kable()` for HMCS table; `moleculaR:::plot.molProbMap()` for top-10 maps

---

## 6 · Dependency Map

```
SIFS
 ├── moleculaR        (spatial windows, SPP, MPM, DSC, FWHM)
 ├── reticulate       (Python bridge for SHAP/XGBoost chunk)
 ├── data.table       (fast CSV I/O for SHAP outputs)
 └── knitr            (vignette rendering, kable tables)

Python (via reticulate)
 ├── xgboost
 ├── shap
 ├── scikit-learn     (GridSearchCV, SelectPercentile)
 ├── joblib           (model serialisation)
 ├── numpy / pandas
 └── matplotlib / seaborn
```

> See `DESCRIPTION` for the formal `Imports` and `Suggests` fields.

---

## 7 · Data Assets (gitignored)

All data assets are hosted on Zenodo: **DOI 10.5281/zenodo.18187395**

| Asset          | Path                                                                          | Description                 |
| -------------- | ----------------------------------------------------------------------------- | --------------------------- |
| Training MIS   | `data/input/0_mis_files/mis_XIII_t3_83318_trainSet.rds`                       | Spatial metadata for t3     |
| Testing MIS    | `data/input/0_mis_files/mis_XV_t2_82558_testSet.rds`                          | Spatial metadata for t4     |
| Training SPMAT | `data/input/01_spmatObjects/trainset_tissue/tr_spmat_XIII_t3_83318.rds`       | Sparse ion matrix t3        |
| Testing SPMAT  | `data/input/01_spmatObjects/testset_tissue/tst_spmat_XV_t1_82558.rds`         | Sparse ion matrix t4        |
| All labels     | `data/input/04_y_lbls/spectraLblsAll/04_Lbl_df_giant_pos/Lbls_df_Pos_all.rds` | Full label dataframe        |
| FWHM object    | `data/input/06_fwhm_object/fwhm_obj_pos.rds`                                  | Instrument peak width model |

Label encoding: **0 = Tumour, 1 = Non-tumour**

---

## 8 · Edit Workflow (Karpathy Pattern)

When the human asks OpenCode to edit any source file:

1. **Read** `wiki/index.md` to locate relevant concept pages.
2. **Read** only the specific pages and source sections that are affected — use `graphify` dependency graph to scope the minimal subgraph.
3. **Discuss** the proposed change with the human before writing.
4. **Edit** the target file using `str_replace` (never rewrite whole files unless necessary).
5. **Update** any impacted wiki concept page.
6. **Append** an entry to `wiki/log.md`: date · file changed · what changed · why.
7. **Update** `wiki/index.md` if new pages were created.

> Use `graphify` to identify the minimal affected subgraph before each edit. Only load the nodes (files/functions) in that subgraph to keep token cost low.

---

## 9 · Vignette Edit Protocol

When the human asks to update `vignettes/01_SIFS_workflow1.Rmd`:

- Treat each numbered section (§5) as an independent edit unit.
- Use `str_replace` on the specific chunk, never overwrite the full file.
- After editing a code chunk, verify that variable names are consistent with the global variables defined in the earlier chunks.
- If a new SIFS function is introduced in the vignette, ensure it exists in `R/` and is documented.
- Check that Python chunks use `r.<varname>` correctly to pass R objects.

---

## 10 · Wiki Page Format

Every wiki page follows this template exactly:

```markdown
# Page Title

**Summary**: One to two sentences.

**Sources**: R source files or vignette chunks this page draws from.

**Last updated**: YYYY-MM-DD

---

Main content. Use clear headings and short paragraphs.

Link related concepts using [[wiki-links]] throughout.

## Related pages

- [[related-concept-1]]
- [[related-concept-2]]
```

---

## 11 · Log & Index Rules

- `wiki/log.md` is **append-only**. Format: `YYYY-MM-DD | <file> | <change summary>`.
- `wiki/index.md` lists every wiki page with a one-line description and last-updated date.
- Both must be updated after every OpenCode operation.

---

## 12 · Rules

| Rule | Detail |
|------|--------|
| Immutable raw | Never modify files in `raw/` |
| GitHub boundary | `wiki/`, `data/`, `raw/`, `output/`, `graphify-out/` are gitignored — never push |
| Minimal edits | Always use `str_replace`; only rewrite whole files if structure changes fundamentally |
| Graphify first | Before any edit, scope the affected subgraph with `graphify` to minimise token use |
| Consistent naming | R objects: `t3_*` / `t4_*` prefix; output files: `<seq>_<descriptor>_<ROI>.<ext>` |
| Log everything | Every operation appended to `wiki/log.md` |
| Language | Plain, precise, scientific — no filler |
| Uncertainty | When ambiguous, ask the human before writing |
