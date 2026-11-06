# Feature extraction and data preparation

The whole path from raw pre-operative DICOM to the PyRadiomics feature matrix: DICOM inventory scan → series-role assignment (T1c identified by `ContrastBolusAgent = Gd-DTPA`) → DICOM→NIfTI conversion → cohort locking and mask-completeness checks → extraction of the 4148 features (Python 3.9 / pyradiomics 3.1.0, fixed bin width 25, no intensity normalisation).

Entry point: `extract_features_pyradiomics.py`

## Files in this folder

- `build_cohort.py`
- `build_vendor_map.py`
- `convert_dicom_to_nifti.py`
- `extract_features_pyradiomics.py`
- `scan_dicom_inventory.py`
- `select_series_roles.py`

### Re-running the necrosis channel (a repair the main batch made necessary)

In the main batch, Raidionics' automatic variant selection for the `MRI_Necrosis` model took the last match by variant-name length. When the input lacked a t1d sequence the selection fell through to the four-sequence `t1c_t1w_t2f_t1d` variant and **silently wrote empty masks**: all 90 Necrosis masks had zero voxels, while the process still returned 0 and the log showed no error. The necrosis channel in the main results tree is therefore unusable and has to be re-run as described below.

- `necrosis_rerun.py` — re-run driver for the single-sequence (t1c) variant. It builds a staging model library `models_nec_t1c/` containing only a `t1c` subdirectory, which forces variant selection onto the single-sequence variant, and writes output to a **separate directory** `patients_nec2/` without touching the validated main results tree. It resumes where it left off: cases that already have a `*Necrosis*.nii.gz` under `T0/` are skipped.
  Usage: `python necrosis_rerun.py [--only <patient>]`
- `necrosis_t1c_pipeline.json` — single-step pipeline definition (model selection only, `MRI_Necrosis`, thresholding reconstruction), read by the script at run time.
- `necrosis_config_template.ini` — Raidionics run-time configuration template (`gpu_id=-1`, `acceleration=cpu`, `reconstruction_order=resample_first`). The placeholders `<WORK_ROOT>` and `<SCRIPT_DIR>` are substituted by the run environment; the script generates the actual configuration per case.
- `run_necrosis_batch.py` — batch launcher: it disables the rename patch inside the re-run script, caps the BLAS thread count, and tees all output to `run_necrosis_batch.log`.
- `necrosis_check.py` — empty-mask checker. It counts label values and non-zero voxels per subregion (TumorCore / Necrosis / Edema / FLAIRChanges), which is how a wholesale channel failure is detected.
  Usage: `python necrosis_check.py <segmentation_tree> [report.txt]`

**Result (n = 90)**: after the re-run, 84/90 cases had a detected necrotic region (median 4.4 mL) and 6 genuinely had none (not a technical failure), against 90/90 empty in the main batch.

> ⚠️ **Data-protection note**: on import, `necrosis_rerun.py` patches `shutil.rmtree`, `os.remove` and `os.rmdir` so that deletions become **renames** (`*_trash_<pid>_<seq>`), which stops Raidionics' internal cleanup from deleting validated data by mistake. Setting `RAIDIONICS_NO_PATCH=1` restores native behaviour, and that is the mode `run_necrosis_batch.py` uses. The patch only renames and never deletes, but it does leave `_trash_*` directories in the output tree, which you can remove yourself afterwards.
