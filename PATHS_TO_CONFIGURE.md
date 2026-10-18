# 需要配置的路径清单 (PATHS_TO_CONFIGURE.md)

本仓库中的脚本保留了对原始分析机器绝对路径的引用，因此**不能直接运行**。
下表列出全部硬编码路径。请把左侧的内容替换为你自己的目录，
推荐做法是在每个脚本顶部把路径常量改为你自己的值，或统一放到一个 `config.py` 中 import。

| 应指向的内容 | 原始路径（示例） | 出现次数 | 涉及脚本 |
|---|---|---|---|
| 原始术前 DICOM 根目录 | `E:\NeuroImagepre` | 4 | `07_acquisition_geometry_qc/audit_vendor_fieldstrength.py`, `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/scan_pre_dicom.py`, `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| 原始术前 DICOM 根目录 | `D:\1桌面文件\AI\NeuroImagepost` | 1 | `01_feature_extraction/build_vendor_map.py` |
| 原始术前 DICOM 根目录 | `D:\NeuroImagepre 全部患者的 DICOM 序列` | 1 | `01_feature_extraction/scan_dicom_inventory.py` |
| 原始术前 DICOM 根目录 | `D:\NeuroImagepre` | 1 | `01_feature_extraction/scan_dicom_inventory.py` |
| 原始术前 DICOM 根目录 | `PATIENT_ROOT`（默认 `D:\NeuroImagepre`）下的 `to_name('P032')\20231228\ImageFileName_*\8` | 1 | `01_feature_extraction/select_series_roles.py` |
| 原始术前 DICOM 根目录 | `PATIENT_ROOT` 下的 `to_name('P051')\20220314\ImageFileName_*\8` | 1 | `01_feature_extraction/select_series_roles.py` |
| 原始术前 DICOM 根目录 | `PATIENT_ROOT` 下的 `to_name('P108')\20230506\ImageFileName_*\12` | 1 | `01_feature_extraction/select_series_roles.py` |
| 原始术前 DICOM 根目录 | `E:\NeuroImagepre\<患者>\ 下的首个可读 DICOM 头，` | 1 | `07_acquisition_geometry_qc/audit_vendor_fieldstrength.py` |
| 原始术前 DICOM 根目录 | `E:\NeuroImagepre（107 例，含 3 例 PAT/STUDY/SERIES 结构）` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| 原始术前 DICOM 根目录 | `E:\\NeuroImagepre` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| 原始术前 DICOM 根目录 | `D:\\NeuroImagepre\\... -> E:\\NeuroImagepre\\...` | 1 | `07_acquisition_geometry_qc/build_s1_paper1_used.py` |
| 原始术前 DICOM 根目录 | `D:\NeuroImagepre\<pat>\<YYYYMMDD>\... 或 D:\NeuroImagepre\<pat>\ImageFileName_xxx\<n>` | 1 | `07_acquisition_geometry_qc/check_study_dates.py` |
| 原始术前 DICOM 根目录 | `E:\NeuroImagepre 为外置盘，本次离线不可读；但既有 S1b（91 例 × 3 序列 = 273 行）` | 1 | `07_acquisition_geometry_qc/rebuild_s1_from_existing.py` |
| 原始术前 DICOM 根目录 | `E:\NeuroImagepre 全量术前 DICOM，提取每个序列的采集参数。` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom.py` |
| 原始术前 DICOM 根目录 | `E:\NeuroImagepre 中 3 个 PAT/STUDY/SERIES/IMG 结构的病例` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| 原始术前 DICOM 根目录 | `to_name("P001")` 所在的 `D:/NeuroImagepre` 子目录 | 1 | `09_figures_and_tables/make_supp_tables_S1_S4.py` |
| Raidionics 分割输出根目录 | `DISC_ROOT`（默认 `D:/Raidionics_data`）下的 `disc_<to_name(Pxxx) 的小写形式>` | 7 | `09_figures_and_tables/make_supp_tables_S1_S4.py` |
| Raidionics 分割输出根目录 | `D:\Raidionics_data\input\2` | 5 | `01_feature_extraction/convert_dicom_to_nifti.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/build_s1_paper1_used.py` … |
| Raidionics 分割输出根目录 | `D:\Raidionics_data\patients\2` | 3 | `01_feature_extraction/build_cohort.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `07_acquisition_geometry_qc/check_s38_claims.py` |
| Raidionics 分割输出根目录 | `D:\Raidionics_data\patients_nec2\2` | 3 | `01_feature_extraction/build_cohort.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `07_acquisition_geometry_qc/check_s38_claims.py` |
| Raidionics 分割输出根目录 | `D:/Raidionics_data/patients/2` | 3 | `06_volume_analysis/make_per_case_volumes.py`, `07_acquisition_geometry_qc/make_s3_slice_thickness.py`, `09_figures_and_tables/make_supp_fig_S1_qc.py` |
| Raidionics 分割输出根目录 | `D:/Raidionics_data/patients_nec2/2` | 3 | `06_volume_analysis/make_per_case_volumes.py`, `07_acquisition_geometry_qc/make_s3_slice_thickness.py`, `09_figures_and_tables/make_supp_fig_S1_qc.py` |
| Raidionics 分割输出根目录 | `D:/Raidionics_data/input/2` | 3 | `07_acquisition_geometry_qc/make_s3_slice_thickness.py`, `09_figures_and_tables/make_supp_fig_S1_qc.py`, `09_figures_and_tables/make_supp_tables_S1_S4.py` |
| Raidionics 分割输出根目录 | `D:\Raidionics_data` | 2 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Raidionics 分割输出根目录 | `D:\Raidionics_data\input\2\<Patient>\T0` | 1 | `01_feature_extraction/convert_dicom_to_nifti.py` |
| Raidionics 分割输出根目录 | `D:\Raidionics_data\disc_*（2 例：P004 / P025）` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Raidionics 分割输出根目录 | `D:\\Raidionics_data` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| Raidionics 分割输出根目录 | `D:\Raidionics_data\disc_* 的 2 例（<PATIENT> / <PATIENT>）。` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| 临床数据表 (GBM_clinical_data.xlsx) | `D:\1桌面文件\AI\GBM_clinical_data.xlsx` | 3 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/build_s1_paper1_used.py`, `07_acquisition_geometry_qc/check_study_vs_surgery.py` |
| 临床数据表 (GBM_clinical_data.xlsx) | `D:/1桌面文件/AI/GBM_clinical_data.xlsx` | 3 | `08_sensitivity_and_power/extended_clinical_baseline.py`, `09_figures_and_tables/make_supp_tables_S1_S4.py`, `09_figures_and_tables/make_table1_baseline.py` |
| 其他临床/合并表 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\series_selection.csv` | 3 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/check_study_dates.py`, `07_acquisition_geometry_qc/diag_series_paths.py` |
| 其他临床/合并表 | `D:\1桌面文件\AI\Glioma_Patients_Classified_merged.xlsx` | 1 | `01_feature_extraction/build_cohort.py` |
| 其他临床/合并表 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\.workbuddy\cohort_res.txt` | 1 | `01_feature_extraction/build_cohort.py` |
| 其他临床/合并表 | `C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts/cohort_MGMT.csv` | 1 | `07_acquisition_geometry_qc/make_s3_slice_thickness.py` |
| 中间结果目录（本仓库的 analysis/） | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1` | 16 | `02_nested_cross_validation/paper1_nested_cv.py`, `03_vendor_robustness/ge_subgroup_permutation_null.py`, `04_deconfounding_decomposition/deconfound_vendor_vs_protocol.py`, `06_volume_analysis/make_per_case_volumes.py` … |
| 中间结果目录（本仓库的 analysis/） | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/analysis` | 3 | `03_vendor_robustness/ge_subgroup_fold_diagnostics.py`, `09_figures_and_tables/make_supp_figs_S3_S4.py`, `09_figures_and_tables/make_supp_tables_S5_S7.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\tables` | 2 | `07_acquisition_geometry_qc/build_acquisition_table.py`, `07_acquisition_geometry_qc/make_s2_radiomics_config.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1` | 2 | `07_acquisition_geometry_qc/build_s1_paper1_used.py`, `07_acquisition_geometry_qc/rebuild_s1_from_existing.py` |
| 中间结果目录（本仓库的 analysis/） | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/figures` | 2 | `09_figures_and_tables/make_supp_fig_S1_qc.py`, `09_figures_and_tables/make_supp_figs_S3_S4.py` |
| 中间结果目录（本仓库的 analysis/） | `D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/tables` | 2 | `09_figures_and_tables/make_supp_tables_S1_S4.py`, `09_figures_and_tables/make_supp_tables_S5_S7.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts` | 1 | `07_acquisition_geometry_qc/build_acquisition_table.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\check_study_date_consistency.csv` | 1 | `07_acquisition_geometry_qc/check_study_dates.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\S1b_full_detailed.csv` | 1 | `07_acquisition_geometry_qc/check_study_vs_surgery.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\check_study_vs_surgery.csv` | 1 | `07_acquisition_geometry_qc/check_study_vs_surgery.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\pre_dicom_series_inventory.csv` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom.py` |
| 中间结果目录（本仓库的 analysis/） | `D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\pre_dicom_series_inventory_extra.csv` | 1 | `07_acquisition_geometry_qc/scan_pre_dicom_extra.py` |
| Python 解释器 / 虚拟环境 | `D:\Anaconda\envs\nnunet 环境运行` | 1 | `01_feature_extraction/extract_features_pyradiomics.py` |
| Python 解释器 / 虚拟环境 | `D:\Anaconda base 运行` | 1 | `05_combat_harmonization/legacy_simplified_combat_exploratory.py` |
| DICOM 或输入暂存目录 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\input_reg\2` | 1 | `01_feature_extraction/extract_features_pyradiomics.py` |
| 本机用户名 | `C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts` | 8 | `03_vendor_robustness/ge_subgroup_fold_diagnostics.py`, `03_vendor_robustness/ge_subgroup_permutation_null.py`, `04_deconfounding_decomposition/deconfound_vendor_vs_protocol.py`, `08_sensitivity_and_power/extended_clinical_baseline.py` … |
| 本机用户名 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts` | 5 | `01_feature_extraction/build_cohort.py`, `01_feature_extraction/extract_features_pyradiomics.py`, `05_combat_harmonization/legacy_simplified_combat_exploratory.py`, `07_acquisition_geometry_qc/build_s1_paper1_used.py` … |
| 本机用户名 | `C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57` | 4 | `02_nested_cross_validation/paper1_nested_cv.py`, `06_volume_analysis/volume_biomarker_analysis.py`, `09_figures_and_tables/make_main_figures.py`, `09_figures_and_tables/make_table1_baseline.py` |
| 本机用户名 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57` | 1 | `01_feature_extraction/build_vendor_map.py` |
| 本机用户名 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\analysis` | 1 | `05_combat_harmonization/legacy_simplified_combat_exploratory.py` |
| 特征目录 / 输出目录（可用环境变量覆盖） | `PAPER1_FEATURE_DIR` / `PAPER1_OUT_DIR` | 1 | `05_combat_harmonization/combat_protected_covariate_sensitivity.py`（默认指向本机 scripts 与 reanalysis 目录；`PAPER1_FEATURE_DIR` 需含 `radiomics_features_long.csv` 与 `cohort_MGMT.csv`） |
| 本机用户名 | `C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\extract_features.py` | 1 | `07_acquisition_geometry_qc/make_s2_radiomics_config.py` |
| 术前影像盘根目录（环境变量可覆盖） | `PATIENT_ROOT`（默认 `D:\NeuroImagepre`） | 1 | `01_feature_extraction/select_series_roles.py` |
| 患者编号 ↔ 本地姓名对照（**含 PHI，不要提交**） | `PAPER1_PATIENT_MAP` 或 `github/tools/patient_ids_local.csv` | 1 | `tools/patient_ids.py`（被 01/02/03/06/07/08/09 共 18 个脚本导入） |
| 临床表拼写变体表（可选） | `PAPER1_PATIENT_ALIAS` 或 `github/tools/patient_alias_local.csv` | 1 | `tools/patient_ids.py` |

---

## 关于患者标识

脚本中**不含任何患者姓名**，患者一律以去标识编号 `P001`–`P109` 表示
（与 Zenodo 数据集、期刊补充表一致）。

- **队列排除**（单例 United Imaging 1.5 T）一律按 `vendor` 判断，不按姓名。
  按姓名排除会在发布版数据上静默失效并把 n 由 89 变成 90。
- **需要在本地运行**时（DICOM/NIfTI 目录以姓名命名），把编号 ↔ 姓名对照放到
  `github/tools/patient_ids_local.csv`（已列入 `.gitignore`），或设置环境变量
  `PAPER1_PATIENT_MAP` 指向该文件。该文件**不得**上传到 GitHub 或 Zenodo。
- 临床表 `GBM_clinical_data.xlsx` 的 `Name` 列若与队列目录名拼写不一致，
  用 `github/tools/patient_alias_local.csv` 补充（模板见同名 `.example` 文件）。
