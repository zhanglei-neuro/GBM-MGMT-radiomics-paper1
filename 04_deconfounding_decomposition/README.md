# Deconfounding decomposition: vendor effect versus acquisition-protocol effect

For every radiomic feature the script fits `z(feature) ~ vendor + 2D/3D + slice thickness + TR + TE + TI + in-plane pixel spacing`, applies BH-FDR to the vendor term, and compares the result pairwise with the vendor-only Kruskal–Wallis test. It also reports partial R² in both directions and the partial R² of in-plane spacing itself.

**Why in-plane pixel spacing is included**: in-plane spacing in this cohort ranges from 0.34 to 1.30 mm and is strongly correlated with vendor (r = 0.375). Omitting the covariate, as an earlier version did, overstates the net vendor effect; with it included, the number of FDR-significant features rises from 259 (6.2%) to 362 (8.7%), the vendor partial R² falls from 0.0164 to 0.0117, and the partial R² of in-plane spacing itself is only 0.0081, most of its information having already been absorbed by the existing acquisition covariates. The conclusion is unchanged.

> Which number is quoted: the manuscript, Table S4 and Figure S4 use the **five-covariate registered specification** (vendor + 2D/3D + slice thickness + TR + TE + TI), which gives 259 (6.2%). The six-covariate model above is a **sensitivity analysis** and is written to the `sensitivity_inplane_covariate` block of the output.

**Where the covariates come from**: the per-case DICOM-measured parameters of the sequence the feature actually belongs to (Supplementary Table `S1b_paper1_per_case.csv`). In-plane spacing is the mean of the row and column spacing.

**Cohort**: the single United Imaging case (n = 1, no platform-level inference possible) is excluded **by vendor**, which keeps n = 89 in line with the manuscript. Excluding by name would break once names are pseudonymised.

Entry point: `deconfound_vendor_vs_protocol.py`

## Files in this folder

- `deconfound_vendor_vs_protocol.py`
