# Acquisition geometry and quality control

Per-case, per-series DICOM acquisition parameter tables (Tables S1a–S1c, S2a–S2d, S3); checks on series role and vendor attribution; verification that the acquisition date precedes the date of surgery; and segmentation-derived metrics stratified by slice thickness.

## Files in this folder

- `audit_vendor_fieldstrength.py`
- `build_acquisition_table.py` — whole-database audit (109 cases, all series). **Its output is not a
  supplementary table**: it includes the single UIH 1.5 T case and therefore does not match the
  89-case cohort. Kept for traceability only; see `supplementary_tables/TABLE_LEGENDS.md`.
- `build_s1_paper1_used.py` — **produces the submitted Tables S1a–S1c** (89 cases, one series per
  sequence; writes `S1a_paper1_by_vendor.csv`, `S1b_paper1_per_case.csv`,
  `S1c_paper1_vendor_signature.csv`)
- `check_s38_claims.py`
- `check_study_dates.py`
- `check_study_vs_surgery.py`
- `diag_series_paths.py`
- `make_s2_radiomics_config.py`
- `make_s3_slice_thickness.py`
- `rebuild_s1_from_existing.py`
- `scan_pre_dicom.py`
- `scan_pre_dicom_extra.py`

> Two scripts write files with the same S1a/S1b/S1c prefixes from different cohorts. Only the
> `*_paper1_*` files from `build_s1_paper1_used.py` are submitted as Tables S1a–S1c; the
> `build_acquisition_table.py` output lives in the private `_internal/` archive.
