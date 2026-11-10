# Figures and tables

Generates the main figures (six in the submitted set: one participant-flow diagram, four single panels and one composite), ESM Figure S5, the main table (Table 1) and the supplementary figures and tables (Figures S1–S4, Tables S1–S7) from the intermediate JSON/CSV results.

## Files in this folder

- `make_figure1_flow.py` — Figure 1, the STARD-AI participant flow diagram (no data input; the
  exclusion counts are hard-coded from Appendix B, section B2)
- `make_main_figures.py`
- `make_composite_main_figures.py` — Figure 3 (ROC + calibration + DCA composite; see
  `figures/FIGURE_LEGENDS.md`)
- `make_figure6_volume_biomarkers.py` — Figure 6, all three panels, rebuilt from the deposited
  Zenodo intermediates (`02_predictions_and_volumes/`, `03_intermediate_results/table3_main.csv`)
  without any raw-DICOM path
- `make_figS5_shap_roi.py` — ESM Figure S5 (SHAP importance plus AUC by ROI × sequence), also
  rebuilt from the deposited intermediates only
- `make_supp_fig_S1_qc.py`
- `make_supp_figs_S3_S4.py`
- `make_supp_tables_S1_S4.py` — writes **Tables S2a–S2d** and the combined workbook
  `Supplementary_Tables_S1_S2.xlsx`. Its S1a/S1b/S1c blocks are the **first-generation** whole-database
  versions and have been superseded by `07_acquisition_geometry_qc/build_s1_paper1_used.py`, which
  produces the submitted Tables S1a–S1c for the 89-case cohort (see
  `supplementary_tables/TABLE_LEGENDS.md`)
- `make_supp_tables_S5_S7.py`
- `make_table1_baseline.py`
