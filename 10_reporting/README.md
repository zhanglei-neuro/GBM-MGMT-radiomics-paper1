# Reporting (submission document generation)

Markdown to Word rendering for the manuscript and the supplementary appendices,
plus the revision log for the text-level corrections applied to the draft.

## Files

- `render_manuscript_docx.py` — generic Markdown to .docx renderer
  (A4, Times New Roman 12 pt, 1.5 line spacing, three-line tables).
- `REVISION_LOG.md` — what the 2026-10-11 text revision changed and why
  (IDH correction, removal of TERT and follow-up content, reference renumbering).

## Usage

```
python render_manuscript_docx.py ../../manuscript/Paper1_Manuscript_MGMT_multivendor_radiomics_EN.md manuscript.docx
```

The renderer understands ATX headings, block quotes, pipe tables, bullet lists,
numbered reference lists and `![caption](file.png)` figure blocks; anything else
is emitted as a body paragraph. Pass `--figdir` when figures are referenced by a
relative path.
