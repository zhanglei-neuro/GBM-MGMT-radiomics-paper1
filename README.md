# Paper 1 — Analysis Code

**Multivendor MRI radiomics for predicting MGMT promoter methylation in glioblastoma: a fully automated segmentation pipeline with explicit vendor-robustness assessment**

This repository contains the complete analysis code for the above manuscript. It covers the whole
pipeline from DICOM inventory to the final figures and tables.

> **Companion data deposit (Zenodo):** https://doi.org/10.5281/zenodo.XXXXXXX
> Per-case radiomic feature matrices, out-of-fold predictions, per-case volumes and all
> intermediate results (JSON/CSV) are archived there. **The data are not in this repository.**

---

## 1. Repository layout

| Folder | Content |
|---|---|
| `01_feature_extraction/` | DICOM inventory, series-role assignment, DICOM→NIfTI conversion, cohort locking, PyRadiomics feature extraction |
| `02_nested_cross_validation/` | **Main analysis**: repeated stratified nested CV (5-fold × 5 repeats), inner feature screening, model comparison (DeLong / NRI / IDI), leave-one-vendor-out (LOVO) |
| `03_vendor_robustness/` | GE-subgroup fold-level diagnostics and label-permutation null distribution |
| `04_deconfounding_decomposition/` | Vendor effect vs. acquisition-protocol effect: per-feature linear model + BH-FDR + partial R² |
| `05_combat_harmonization/` | ComBat vendor harmonization — **the implementation reported in the paper is inside `02_nested_cross_validation/paper1_nested_cv.py` §[7] (neuroCombat 0.2.12, empirical-Bayes shrinkage)**; this folder only keeps a superseded exploratory script |
| `06_volume_analysis/` | Seven pre-specified volume biomarkers (native-grid voxel counting) and their association with MGMT status |
| `07_acquisition_geometry_qc/` | Per-case/per-series acquisition parameter tables (Tables S1a–S1c, S2a–S2d, S3) and geometry consistency checks |
| `08_sensitivity_and_power/` | MDES power analysis (Hanley–McNeil), stronger clinical baselines, cohort- and geometry-stratified sensitivity analyses |
| `09_figures_and_tables/` | Generation of main Figures 1–8, main Table 1, supplementary Figures S1–S4 and Tables S1–S7 |
| `10_reporting/` | Markdown → Word rendering of the manuscript and the response letter; revision patches |
| `tools/` | Safety and reproducibility helpers (see §4) |
| `tools/patient_ids.py` | Patient-identifier normalisation — **required by most scripts**, see §4 |

## 2. Script index

### 01_feature_extraction
| Script | Purpose |
|---|---|
| `scan_dicom_inventory.py` | Scan the pre-operative DICOM tree and build a series inventory |
| `select_series_roles.py` | Assign each series a role (T1c / T1w / FLAIR) using DICOM header rules (`ContrastBolusAgent = Gd-DTPA` for T1c; `dark-fluid`/`tirm` for FLAIR) |
| `convert_dicom_to_nifti.py` | Rebuild the Raidionics input directory and convert DICOM → NIfTI (role-suffixed filenames) |
| `build_cohort.py` | Lock the modelling cohort and verify mask completeness (TumorCore / Edema / FLAIRChanges / Necrosis) |
| `build_vendor_map.py` | Build the per-case scanner-vendor map from DICOM headers |
| `extract_features_pyradiomics.py` | **PyRadiomics feature extraction** (Python 3.9 / pyradiomics 3.1.0): 4 ROIs × 3 sequences = 4148 features per case, fixed bin width 25, no intensity normalization |

### 02_nested_cross_validation
| Script | Purpose |
|---|---|
| `paper1_nested_cv.py` | Main analysis. Global label-agnostic pre-processing (missing >20 % and zero variance), repeated stratified nested CV (outer 5-fold × 5 repeats; inner: univariate \|AUC−0.5\| top-200 → greedy de-correlation \|Spearman\|>0.9 → ≤40 features → 3-fold hyper-parameter search). Models: L2 logistic regression, XGBoost, random forest; controls: clinical baseline (age + sex) and combined model. Reports repeat-mean ± SD and pooled-OOF AUC, DeLong test, category-free NRI, IDI, calibration slope/intercept, Brier score, bootstrap 95 % CI, and LOVO. All pre-processing is fitted inside the training folds only. |

### 03_vendor_robustness
| Script | Purpose |
|---|---|
| `ge_subgroup_fold_diagnostics.py` | Per-fold class distribution, per-fold AUC and direction stability inside the GE subgroup (n = 32) |
| `ge_subgroup_permutation_null.py` | 100 label permutations with the identical pipeline, to test whether the sub-chance GE AUC is explicable by chance |

### 04_deconfounding_decomposition
| Script | Purpose |
|---|---|
| `deconfound_vendor_vs_protocol.py` | Fits `z(feature) ~ vendor + 2D/3D + slice thickness + TR + TE + TI + in-plane pixel spacing (mean of row/col)` for each of the 4148 features; BH-FDR on the vendor term; paired comparison with the vendor-only Kruskal–Wallis test; partial R² in both directions and for the in-plane spacing term. Cohort excludes the single United Imaging case by vendor (n = 89). |

### 05_combat_harmonization
| Script | Purpose |
|---|---|
| *(the reported method lives in `02_nested_cross_validation/paper1_nested_cv.py` §[7])* | ComBat harmonization as reported in the paper: **neuroCombat 0.2.12** (empirical-Bayes shrinkage), fitted inside each training fold and applied to the test fold via `neuroCombatFromTraining` (no leakage); batch = vendor only, no biological covariates; plus the transductive-free comparator VWS (vendor-wise Z-standardisation with training-fold-only parameters) |
| `combat_protected_covariate_sensitivity.py` | **Reported in §3.7**: sensitivity analysis repeating the same nested CV with MGMT added as a protected discrete covariate (`categorical_cols=['MGMT']`) inside the training folds. Protecting the outcome-associated variance did **not** restore discrimination (logistic ΔAUC −0.017, p = 0.77; XGBoost ΔAUC −0.127, p = 0.09) |
| `legacy_simplified_combat_exploratory.py` | **Superseded / not reported.** Simplified location/scale ComBat (Johnson 2007, no empirical-Bayes shrinkage) as a fit/transform step. Kept for provenance only; its outputs are not used in the paper |

### 06_volume_analysis
| Script | Purpose |
|---|---|
| `make_per_case_volumes.py` | Per-case volumes (TC / ED / FL / BR / NE) by **native-grid voxel counting**, avoiding resampling drift |
| `volume_biomarker_analysis.py` | The seven pre-specified volume biomarkers, Mann–Whitney tests, univariate AUC with bootstrap CI, within-family BH-FDR, and the multivariable model |

### 07_acquisition_geometry_qc
| Script | Purpose |
|---|---|
| `build_acquisition_table.py` | Per-series acquisition parameter table (Tables S1a–S1c) |
| `audit_vendor_fieldstrength.py` | Vendor and field-strength audit |
| `make_s2_radiomics_config.py` | PyRadiomics configuration and the derivation of the 4148-dimensional feature space (Tables S2a–S2d) |
| `make_s3_slice_thickness.py` | Segmentation-derived metrics stratified by T1c acquisition geometry (Table S3) |
| `build_s1_paper1_used.py`, `rebuild_s1_from_existing.py` | Restrict S1 to the series actually used for radiomics |
| `scan_pre_dicom.py`, `scan_pre_dicom_extra.py` | Pre-operative DICOM scan inventory (including PAT/STUDY/SERIES structures) |
| `check_study_dates.py`, `check_study_vs_surgery.py` | Verify imaging date < surgery date for every case |
| `diag_series_paths.py`, `check_s38_claims.py` | Series-path diagnostics; verification of the claims made in §3.8 |

### 08_sensitivity_and_power
| Script | Purpose |
|---|---|
| `power_mdes_and_review_analyses.py` | **Minimum detectable AUC (MDES)** under α = 0.05 two-sided, power = 0.80/0.90 (Hanley–McNeil; bi-normal equal-variance model), plus the analyses added during peer review |
| `extended_clinical_baseline.py` | Stronger clinical baselines (age + sex + temporal-lobe involvement + tumour volume ± oedema volume) and calibration slope/intercept |
| `sensitivity_single_repeat.py` | Sensitivity cohorts A–D using a single-repeat pooled-OOF protocol |
| `sensitivity_exclude_uih_qc.py` | Exclusion of the United Imaging 1.5 T single case; stratification by acquisition geometry |
| `merge_sensitivity_results.py` | Merge sensitivity results into the supplementary tables |

### 09_figures_and_tables
| Script | Output |
|---|---|
| `make_main_figures.py` | Figures 1–8 |
| `make_table1_baseline.py` | Table 1 (baseline, by vendor) |
| `make_supp_tables_S5_S7.py` | Tables S5–S7 |
| `make_supp_tables_S1_S4.py` | Tables S1–S4 |
| `make_supp_figs_S3_S4.py` | Figures S3–S4 |
| `make_supp_fig_S1_qc.py` | Figure S1 (segmentation QC) |

### 10_reporting
| Script | Purpose |
|---|---|
| `render_manuscript_docx.py`, `render_response_docx.py` | Markdown → Word rendering of the manuscript and the response letter |
| `patch_manuscript_revision.py`, `patch_response_revision.py` | Textual revision patches applied during peer review |

## 3. Environment

Two environments are required.

**A. Python 3.12 — statistics, modelling, figures, tables**

```bash
pip install -r requirements.txt
```

Key versions used for the reported results:
`scikit-learn 1.5.1`, `xgboost 3.4.1`, `numpy 1.26.4`, `pandas 2.x`, `shap`,
`neuroCombat 0.2.12`, `matplotlib`, `statsmodels`.

**B. Python 3.9 — PyRadiomics feature extraction only**

```bash
conda create -n nnunet python=3.9
conda activate nnunet
pip install pyradiomics==3.1.0 pydicom SimpleITK numpy
```

Random seed is fixed at **42** throughout.

## 4. Patient identifier convention (`tools/patient_ids.py`)

The released data use the sequential study code **`P001`–`P109`**. The scripts follow three rules:

1. **Data keys are always `Pxxx`.** After loading any patient-indexed table the identifier is
   normalised with `to_id()`, so the same code works on the released data (codes) and on a local
   working copy (names).
2. **Filesystem paths use the local directory names.** DICOM/NIfTI directories are named after the
   patient; paths therefore go through `to_name(pid)`.
3. **Cohort exclusions are decided by `vendor`, never by name.** The single United Imaging 1.5 T
   case is removed with `vendor != 'UIH'`. Excluding by name would silently fail on the released
   data (where names are codes) and change *n* from 89 to 90.

For a local run, provide the code ↔ name mapping (kept outside the repository):

```bash
cp <your secure location>/patient_ids_local.csv   github/tools/patient_ids_local.csv
# or:  export PAPER1_PATIENT_MAP=/secure/path/patient_ids_local.csv
```

`tools/patient_ids_local.csv` is git-ignored. Spelling variants between the clinical spreadsheet and
the cohort directory names are supplied by the optional, PHI-free
`tools/patient_alias_local.csv` (template: `tools/patient_alias_local.csv.example`).

## 4b. Tools

| Tool | Purpose |
|---|---|
| `tools/check_phi.py` | **Run this before every push.** Scans the repository for patient identifiers (names, Chinese names, original DICOM paths) and fails loudly if any are found |
| `tools/patient_ids.py` | Identifier normalisation + vendor-based cohort exclusion (see §4) |
| `tools/sanitize_phi.py` | Replaces any residual patient identifier with the `<PATIENT>` placeholder |
| `tools/configure_paths.py` | Interactively rewrites the hard-coded local paths in every script to your own environment |

> ⚠️ **Patient privacy.** The repository contains **no** patient names: every patient-level
> reference is the sequential code `P001`–`P109` used in the archived data, and the name ↔ code
> mapping (`tools/patient_ids_local.csv`, `tools/patient_alias_local.csv`) is deliberately **not**
> part of the repository. Never commit DICOM files, the clinical spreadsheet, or any file
> containing patient names.

## 5. Path configuration

The scripts were written against the original analysis machine and contain absolute paths
(Windows drive letters). **They therefore do not run out of the box.** Run

```bash
python tools/configure_paths.py
```

or edit the constants at the top of each script. `PATHS_TO_CONFIGURE.md` lists every
hard-coded path and what it should point to.

## 6. Expected inputs / outputs

| Stage | Input | Output |
|---|---|---|
| 01 | pre-operative DICOM tree | `radiomics_features_long.csv` (Zenodo `01_feature_matrix/`) |
| 02 | feature matrix + labels | `predictions.csv`, `metrics_main.json` (Zenodo `02_…`, `03_…`) |
| 03 | feature matrix + labels + vendor | `ge_diagnostics.json`, `ge_perm_null.json` |
| 04 | feature matrix + per-series DICOM parameters | `m4_vendor_vs_protocol.json` |
| 05 | feature matrix + batch labels | `metrics_combat.json` |
| 06 | segmentation masks | `per_case_volumes.csv`, `metrics_volume.json` |
| 08 | feature matrix + labels | `review_metrics.json`, `sens_cohort_*.json` |
| 09 | all of the above | Figures 1–8, Tables 1–3, S1–S7 |

## 7. Citation

```bibtex
@article{paper1_mgmt_radiomics,
  title   = {Multivendor MRI radiomics for predicting MGMT promoter methylation in
             glioblastoma: a fully automated segmentation pipeline with explicit
             vendor-robustness assessment},
  author  = {<Author list>},
  journal = {<Journal>},
  year    = {2026},
  note    = {Code: https://github.com/<user>/<repo>; Data: https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

## 8. Licence

MIT — see `LICENSE`. The archived data on Zenodo are released under CC BY 4.0.
