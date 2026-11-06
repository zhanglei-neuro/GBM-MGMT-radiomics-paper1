# Paths that must be configured (PATHS_TO_CONFIGURE.md)

The scripts keep references to absolute paths on the original analysis machine, so **they do not run as they stand**. The table below lists every hard-coded path. Replace what is in the left-hand column with your own directory — preferably by editing the path constants at the top of each script, or by putting them all in a `config.py` that the scripts import.

> A few example paths below contain a Chinese folder name, `桌面文件` (literally "desktop files"), because that is the real directory name on the original workstation. They are quoted verbatim so that you can still find them by searching the code; you will be replacing them in any case.

| What it should point to | Original path (example) | Occurrences | Scripts |
|---|---|---|---|
| Root of the original pre-operative DICOM | `E:\NeuroImagepre` | 4 | `07_acquisition_geometry_qc/audit_vendor_fieldstrength.py`, `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/scan_pre_dicom.py`, `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Root of the original pre-operative DICOM | `D:\1桌面文件\AI\NeuroImagepost` | 1 | `01_feature_extraction/build_vendor_map.py` |
| Root of the original pre-operative DICOM | `D:\NeuroImagepre` (the DICOM tree of all patients) | 1 | `01_feature_extraction/scan_dicom_inventory.py` |
| Root of the original pre-operative DICOM | `D:\NeuroImagepre` | 1 | `01_feature_extraction/scan_dicom_inventory.py` |
| Root of the original pre-operative DICOM | `to_name('P032')\20231228\ImageFileName_*\8` under `PATIENT_ROOT` (default `D:\NeuroImagepre`) | 1 | `01_feature_extraction/select_series_roles.py` |
| Root of the original pre-operative DICOM | `to_name('P051')\20220314\ImageFileName_*\8` under `PATIENT_ROOT` | 1 | `01_feature_extraction/select_series_roles.py` |
| Root of the original pre-operative DICOM | `to_name('P108')\20230506\ImageFileName_*\12` under `PATIENT_ROOT` | 1 | `01_feature_extraction/select_series_roles.py` |
| Root of the original pre-operative DICOM | the first readable DICOM header under `E:\NeuroImagepre\<patient>\` | 1 | `07_acquisition_geometry_qc/audit_vendor_fieldstrength.py` |
| Root of the original pre-operative DICOM | `E:\NeuroImagepre` (107 cases, including 3 with a PAT/STUDY/SERIES structure) | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Root of the original pre-operative DICOM | `E:\\NeuroImagepre` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Root of the original pre-operative DICOM | `D:\\NeuroImagepre\\... -> E:\\NeuroImagepre\\...` | 1 | `07_acquisition_geometry_qc/build_s1_paper1_used.py` |
| Root of the original pre-operative DICOM | `D:\NeuroImagepre\<pat>\<YYYYMMDD>\...` or `D:\NeuroImagepre\<pat>\ImageFileName_xxx\<n>` | 1 | `07_acquisition_geometry_qc/check_study_dates.py` |
| Root of the original pre-operative DICOM | `E:\NeuroImagepre` (external drive, unreadable offline at that point; the existing S1b has 91 cases × 3 sequences = 273 rows) | 1 | `07_acquisition_geometry_qc/rebuild_s1_from_existing.py` |
| Root of the original pre-operative DICOM | the full pre-operative DICOM tree under `E:\NeuroImagepre`, from which acquisition parameters are extracted for each series | 1 | `07_acquisition_geometry_qc/scan_pre_dicom.py` |
| Root of the original pre-operative DICOM | the 3 cases under `E:\NeuroImagepre` that use a PAT/STUDY/SERIES/IMG structure | 1 | `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Root of the original pre-operative DICOM | the `D:/NeuroImagepre` subdirectory holding `to_name("P001")` | 1 | `09_figures_and_tables/make_supp_tables_S1_S4.py` |
| Raidionics segmentation output root | `disc_<lower-case to_name(Pxxx)>` under `DISC_ROOT` (default `D:/Raidionics_data`) | 7 | `09_figures_and_tables/make_supp_tables_S1_S4.py` |
| Raidionics segmentation output root | `D:\Raidionics_data\input\2` | 5 | `01_feature_extraction/convert_dicom_to_nifti.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/build_s1_paper1_used.py`, … |
| Raidionics segmentation output root | `D:\Raidionics_data\patients\2` | 3 | `01_feature_extraction/build_cohort.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `07_acquisition_geometry_qc/check_s38_claims.py` |
| Raidionics segmentation output root | `D:\Raidionics_data\patients_nec2\2` | 3 | `01_feature_extraction/build_cohort.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `07_acquisition_geometry_qc/check_s38_claims.py` |
| Raidionics segmentation output root | `D:/Raidionics_data/patients/2` | 3 | `06_volume_analysis/make_per_case_volumes.py`, `07_acquisition_geometry_qc/make_s3_slice_thickness.py`, `09_figures_and_tables/make_supp_fig_S1_qc.py` |
| Raidionics segmentation output root | `D:/Raidionics_data/patients_nec2/2` | 3 | `06_volume_analysis/make_per_case_volumes.py`, `07_acquisition_geometry_qc/make_s3_slice_thickness.py`, `09_figures_and_tables/make_supp_fig_S1_qc.py` |
| Raidionics segmentation output root | `D:/Raidionics_data/input/2` | 3 | `07_acquisition_geometry_qc/make_s3_slice_thickness.py`, `09_figures_and_tables/make_supp_fig_S1_qc.py`, `09_figures_and_tables/make_supp_tables_S1_S4.py` |
| Raidionics segmentation output root | `D:\Raidionics_data` | 2 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Raidionics segmentation output root | `D:\Raidionics_data\input\2\<Patient>\T0` | 1 | `01_feature_extraction/convert_dicom_to_nifti.py` |
| Raidionics segmentation output root | `D:\Raidionics_data\disc_*` (2 cases: P004 / P025) | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Raidionics segmentation output root | `D:\\Raidionics_data` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Raidionics segmentation output root | the 2 cases under `D:\Raidionics_data\disc_*` (`<PATIENT>` / `<PATIENT>`) | 1 | `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Clinical data table (GBM_clinical_data.xlsx) | `D:\1桌面文件\AI\GBM_clinical_data.xlsx` | 3 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/build_s1_paper1_used.py`, `07_acquisition_geometry_qc/check_study_vs_surgery.py` |
| Clinical data table (GBM_clinical_data.xlsx) | `D:/1桌面文件/AI/GBM_clinical_data.xlsx` | 3 | `08_sensitivity_and_power/extended_clinical_baseline.py`, `09_figures_and_tables/make_supp_tables_S1_S4.py`, `09_figures_and_tables/make_table1_baseline.py` |
| Other clinical and merged tables | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\series_selection.csv` | 3 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/check_study_dates.py`, `07_acquisition_geometry_qc/diag_series_paths.py` |
| Other clinical and merged tables | `D:\1桌面文件\AI\Glioma_Patients_Classified_merged.xlsx` | 1 | `01_feature_extraction/build_cohort.py` |
| Other clinical and merged tables | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\.workbuddy\cohort_res.txt` | 1 | `01_feature_extraction/build_cohort.py` |
| Other clinical and merged tables | `C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts/cohort_MGMT.csv` | 1 | `07_acquisition_geometry_qc/make_s3_slice_thickness.py` |
| Intermediate-results directory (the `analysis/` folder of this repository) | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1` | 16 | `02_nested_cross_validation/paper1_nested_cv.py`, `03_vendor_robustness/ge_subgroup_permutation_null.py`, `04_deconfounding_decomposition/deconfound_vendor_vs_protocol.py`, `06_volume_analysis/make_per_case_volumes.py`, … |
| Intermediate-results directory | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/analysis` | 3 | `03_vendor_robustness/ge_subgroup_fold_diagnostics.py`, `09_figures_and_tables/make_supp_figs_S3_S4.py`, `09_figures_and_tables/make_supp_tables_S5_S7.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\tables` | 2 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/make_s2_radiomics_config.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1` | 2 | `07_acquisition_geometry_qc/build_s1_paper1_used.py`, `07_acquisition_geometry_qc/rebuild_s1_from_existing.py` |
| Intermediate-results directory | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/figures` | 2 | `09_figures_and_tables/make_supp_fig_S1_qc.py`, `09_figures_and_tables/make_supp_figs_S3_S4.py` |
| Intermediate-results directory | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/tables` | 2 | `09_figures_and_tables/make_supp_tables_S1_S4.py`, `09_figures_and_tables/make_supp_tables_S5_S7.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\check_study_date_consistency.csv` | 1 | `07_acquisition_geometry_qc/check_study_dates.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\S1b_full_detailed.csv` | 1 | `07_acquisition_geometry_qc/check_study_vs_surgery.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\check_study_vs_surgery.csv` | 1 | `07_acquisition_geometry_qc/check_study_vs_surgery.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\pre_dicom_series_inventory.csv` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom.py` |
| Intermediate-results directory | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\pre_dicom_series_inventory_extra.csv` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Python interpreter / virtual environment | run under the `D:\Anaconda\envs\nnunet` environment | 1 | `01_feature_extraction/extract_features_pyradiomics.py` |
| Python interpreter / virtual environment | run under the `D:\Anaconda` base environment | 1 | `05_combat_harmonization/legacy_simplified_combat_exploratory.py` |
| DICOM or input staging directory | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\input_reg\2` | 1 | `01_feature_extraction/extract_features_pyradiomics.py` |
| Local user name (home directory) | `C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts` | 8 | `03_vendor_robustness/ge_subgroup_fold_diagnostics.py`, `03_vendor_robustness/ge_subgroup_permutation_null.py`, `04_deconfounding_decomposition/deconfound_vendor_vs_protocol.py`, `08_sensitivity_and_power/extended_clinical_baseline.py`, … |
| Local user name (home directory) | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts` | 5 | `01_feature_extraction/build_cohort.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `05_combat_harmonization/legacy_simplified_combat_exploratory.py`, `07_acquisition_geometry_qc/build_s1_paper1_used.py`, … |
| Local user name (home directory) | `C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57` | 4 | `02_nested_cross_validation/paper1_nested_cv.py`, `06_volume_analysis/volume_biomarker_analysis.py`, `09_figures_and_tables/make_main_figures.py`, `09_figures_and_tables/make_table1_baseline.py` |
| Local user name (home directory) | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57` | 1 | `01_feature_extraction/build_vendor_map.py` |
| Local user name (home directory) | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\analysis` | 1 | `05_combat_harmonization/legacy_simplified_combat_exploratory.py` |
| Feature directory / output directory (can be overridden by environment variables) | `PAPER1_FEATURE_DIR` / `PAPER1_OUT_DIR` | 1 | `05_combat_harmonization/combat_protected_covariate_sensitivity.py` (defaults to the local `scripts` and `reanalysis` directories; `PAPER1_FEATURE_DIR` must contain `radiomics_features_long.csv` and `cohort_MGMT.csv`) |
| Local user name (home directory) | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\extract_features.py` | 1 | `07_acquisition_geometry_qc/make_s2_radiomics_config.py` |
| Root of the pre-operative imaging drive (can be overridden by environment variable) | `PATIENT_ROOT` (default `D:\NeuroImagepre`) | 1 | `01_feature_extraction/select_series_roles.py` |
| Patient code ↔ local name mapping (**contains PHI — do not commit**) | `PAPER1_PATIENT_MAP` or `github/tools/patient_ids_local.csv` | 1 | `tools/patient_ids.py` (imported by 18 scripts across 01/02/03/06/07/08/09) |
| Clinical-table spelling variants (optional) | `PAPER1_PATIENT_ALIAS` or `github/tools/patient_alias_local.csv` | 1 | `tools/patient_ids.py` |

---

## About patient identifiers

The scripts contain **no patient names**. Every patient appears as the de-identified code `P001`–`P109`, matching the Zenodo dataset and the journal's supplementary tables.

- **Cohort exclusion** (the single United Imaging 1.5 T case) is always decided by `vendor`, never by name. Excluding by name silently fails on the released data and turns n from 89 into 90.
- **To run the code locally** (where DICOM and NIfTI directories are named after patients), put the code ↔ name mapping at `github/tools/patient_ids_local.csv` (already in `.gitignore`), or set the environment variable `PAPER1_PATIENT_MAP` to point at it. That file **must not** be uploaded to GitHub or Zenodo.
- If the `Name` column of `GBM_clinical_data.xlsx` spells a patient differently from the cohort directory, add it through `github/tools/patient_alias_local.csv` (template: the `.example` file beside it).
