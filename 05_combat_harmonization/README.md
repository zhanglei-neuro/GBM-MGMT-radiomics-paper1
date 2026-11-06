# ComBat harmonisation

> **This folder does not contain the implementation of the ComBat results reported in the paper.** The real implementation lives inside the main analysis pipeline: **§[7] `ComBat harmonisation (neuroCombat)` of `../02_nested_cross_validation/paper1_nested_cv.py`**.

## The method as reported (authoritative definition)

- Implementation: **neuroCombat 0.2.12** (empirical-Bayes shrinkage; Johnson 2007 / Fortin 2017)
- Fitting: inside each outer fold's **training folds**, with the test fold transformed through `neuroCombatFromTraining`, so no parameters leak
- Batch variable: the vendor label only; **no biological covariate is passed in** (the MGMT label plays no part in harmonisation at any point)
- What is reported: nested CV after harmonisation (logistic / XGBoost), LOVO after harmonisation, the number of vendor-significant features after harmonisation against before, and LOVO for the transductive-free alternative **VWS** (vendor-wise Z-standardisation with mean and variance estimated in the training folds only)
- Output: `zenodo/03_intermediate_results/metrics_combat.json`

To reproduce:

```bash
python 02_nested_cross_validation/paper1_nested_cv.py
```

## Files in this folder

- `combat_protected_covariate_sensitivity.py` — **the protected-covariate sensitivity analysis reported in §3.7 of the manuscript.** Within the same nested CV framework it compares two arms: A passes only the vendor batch variable (the main-analysis setting); B additionally passes MGMT as a discrete covariate inside the training folds (`categorical_cols=['MGMT']`). Both are fitted in the training folds and the test fold is transformed with `neuroCombatFromTraining`. Result: protecting the covariate did **not** restore discrimination (logistic ΔAUC −0.017, p = 0.77; XGBoost ΔAUC −0.127, p = 0.09). Output: `zenodo/03_intermediate_results/combat_protected_covariate_sensitivity.json`.
- `legacy_simplified_combat_exploratory.py` — **early exploratory implementation, not used in the manuscript.** It is the simplified location/scale version without empirical-Bayes shrinkage (Johnson 2007), usable as fit/transform, and downstream it reuses the same nested CV pipeline as the main analysis. It is kept for provenance only; its outputs (`predictions_combat.csv` and the older `metrics_combat.json`) were **not used**, and should not be used to reproduce the reported results.

> Note: the file here was originally called `combat_harmonization.py`, which was easily confused with the method the paper reports. It was renamed and given the deprecation banner in the 2026-10-12 revision.
