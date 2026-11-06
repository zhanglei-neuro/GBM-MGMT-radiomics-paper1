# Volume biomarkers

Seven prespecified volumetric markers: the volumes of TumorCore, Necrosis, Edema and FLAIRChanges, total tumour volume, the total tumour to oedema ratio, and the enhancing core to oedema ratio. Volumes come from **voxel counting on the native grid** of each mask, which avoids the volume drift introduced by resampling. Statistics: Mann–Whitney, univariable AUC with bootstrap CI, within-family BH-FDR.

## Files in this folder

- `make_per_case_volumes.py`
- `volume_biomarker_analysis.py`
