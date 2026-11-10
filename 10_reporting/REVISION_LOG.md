# Revision log — text-level corrections applied on 2026-10-11

This folder holds the tooling that turns the Markdown source into submission-ready
Word files. The English package previously also carried three scripts that patched
the **Chinese** draft (`patch_manuscript_revision.py`, `patch_response_revision.py`,
`render_response_docx.py`). Those scripts operate on Chinese anchor sentences and
therefore cannot be carried into an English-only package; they are replaced by this
log, which records exactly what they changed so the revision remains auditable.

The English manuscript in `../../manuscript/` already contains all of these
corrections, so no patch step is needed to reproduce it.

## 1. Scope restriction: no TERT, no follow-up endpoints

Paper 1 is restricted to MGMT promoter methylation. TERT promoter status and all
follow-up endpoints (overall survival, progression-free survival, recurrence) were
removed from the abstract, introduction, methods, results, discussion, conclusion,
Table 1 footnotes and the supplementary material.

Reference [23] of the previous draft dealt with TERT and was deleted; references
24–31 were renumbered to 23–30, and every in-text citation was updated accordingly.

## 2. Acquisition-geometry caveat (Results 3.8)

The thin-slice subgroup restriction was rewritten to state that the 3D thin-slice
share differs between vendors (Siemens 49.1 % vs GE 96.9 %, p < 0.001) and that this
imbalance is the one confirmed vendor-related confounder in the cohort, rather than
a limitation of the segmentation.

## 3. Response letter

The same groups of edits were applied to the point-by-point response so that
the letter and the manuscript quote identical numbers and identical wording.

## 4. Word rendering

`render_manuscript_docx.py` in this folder renders any of the English Markdown files
to .docx (A4, Times New Roman 12 pt, 1.5 line spacing, three-line tables). It
replaces the earlier Chinese-typography renderer, which hard-coded SimSun / SimHei and
the Chinese draft path:

```
python render_manuscript_docx.py ../../manuscript/Paper1_Manuscript_MGMT_multivendor_radiomics_EN.md out.docx
```
