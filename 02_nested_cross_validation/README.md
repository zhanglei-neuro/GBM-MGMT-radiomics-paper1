# Nested cross-validation (main analysis)

The main analysis of the paper. Global label-agnostic pre-processing → outer stratified 5-fold × 5 repeats → inner univariate |AUC−0.5| top-200 → greedy de-correlation at |Spearman| > 0.9 (≤40 features) → 3-fold hyper-parameter search. All pre-processing is fitted in the training folds only.
Models: L2 logistic regression, XGBoost, random forest; comparators: clinical baseline (age + sex) and the combined model.
Metrics: repeat mean ± SD, pooled OOF AUC, DeLong test, category-free NRI, IDI, calibration slope and intercept, Brier score, bootstrap 95% CI, LOVO.

Entry point: `paper1_nested_cv.py`

## Files in this folder

- `paper1_nested_cv.py`
