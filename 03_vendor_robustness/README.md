# Vendor robustness

Fold-by-fold diagnostics inside the GE subgroup (n = 32): class distribution in the training and test folds, per-fold AUC and direction, and the stability of the inner-screening direction; plus a 100-permutation label-shuffle null distribution, used to judge whether the below-chance AUC inside the subgroup can be put down to random fluctuation.

## Files in this folder

- `ge_subgroup_fold_diagnostics.py`
- `ge_subgroup_permutation_null.py`
