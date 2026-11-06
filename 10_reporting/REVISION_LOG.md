# Revision log — text-level corrections applied on 2026-10-11

This folder holds the tooling that turns the Markdown source into submission-ready
Word files. The English package previously also carried three scripts that patched
the **Chinese** draft (`patch_manuscript_revision.py`, `patch_response_revision.py`,
`render_response_docx.py`). Those scripts operate on Chinese anchor sentences and
therefore cannot be carried into an English-only package; they are replaced by this
log, which records exactly what they changed so the revision remains auditable.

The English manuscript in `../../manuscript/` already contains all of these
corrections, so no patch step is needed to reproduce it.

## 1. IDH status (author decision)

The source database had recorded 6 cases as IDH-mutant. This was an entry error.
All 89 cases in the cohort are **IDH-wildtype GBM** as defined by WHO 2021.

Consequences propagated through the text:

| Location | Change |
|---|---|
| Abstract | "89 GBM cases" reworded to state explicitly that the whole cohort is IDH-wildtype GBM (WHO 2021) |
| Introduction | "83 IDH-wildtype GBM + 6 IDH-mutant grade 4 astrocytoma" replaced by a single entity, 89 IDH-wildtype GBM |
| Methods 2.1 (diagnostic criteria) | Cohort definition restated as one entity; no IDH subgroup |
| Methods 2.2 (baseline comparability) | IDH row removed from the comparability table |
| Methods 2.8 (missing data) | IDH missingness paragraph removed |
| Results 3.1, 3.8 | IDH-stratified numbers removed; the sensitivity-analysis row "E" deleted from the 3.8 table |
| Discussion 4.2, 4.4 | IDH-stratification limitation removed; limitation numbering renumbered |
| Conclusion | Sentence on the IDH-mutant subgroup removed |
| Figure/table list | Row referring to the IDH sensitivity analysis deleted |

## 2. Scope restriction: no TERT, no follow-up endpoints

Paper 1 is restricted to MGMT promoter methylation. TERT promoter status and all
follow-up endpoints (overall survival, progression-free survival, recurrence) were
removed from the abstract, introduction, methods, results, discussion, conclusion,
Table 1 footnotes and the supplementary material.

Reference [23] of the previous draft dealt with TERT and was deleted; references
24–31 were renumbered to 23–30, and every in-text citation was updated accordingly.

## 3. Acquisition-geometry caveat (Results 3.8)

The thin-slice subgroup restriction was rewritten to state that the 3D thin-slice
share differs between vendors (Siemens 49.1 % vs GE 96.9 %, p < 0.001) and that this
imbalance is the one confirmed vendor-related confounder in the cohort, rather than
a limitation of the segmentation.

## 4. Response letter

The same four groups of edits were applied to the point-by-point response so that
the letter and the manuscript quote identical numbers and identical wording.

## 5. Word rendering

`render_manuscript_docx.py` in this folder renders any of the English Markdown files
to .docx (A4, Times New Roman 12 pt, 1.5 line spacing, three-line tables). It
replaces the earlier Chinese-typography renderer, which hard-coded SimSun / SimHei and
the Chinese draft path:

```
python render_manuscript_docx.py ../../manuscript/Paper1_Manuscript_MGMT_multivendor_radiomics_EN.md out.docx
```
