# Tools

`patient_ids.py` — normalises patient identifiers (public code ↔ local name) and excludes the single case by vendor;
`check_phi.py` — pre-push privacy self-check; `sanitize_phi.py` — replaces names with the placeholder;
`configure_paths.py` — rewrites the hard-coded paths in bulk.

**Run `check_phi.py` before every `git push`.**

## Files in this folder

- `patient_ids.py` — imported by 18 scripts across 01/02/03/06/07/08/09; see §4 of the repository README
- `patient_alias_local.csv.example` — template for the clinical-table spelling variants (optional)
- `check_phi.py`
- `configure_paths.py`
- `sanitize_phi.py`

## Lookup tables needed for a local run (**do not commit**)

| File | Content | How to obtain it |
|---|---|---|
| `patient_ids_local.csv` | two columns, `patient_id_public,patient_name_original` | copy it here from your secure location, or point `PAPER1_PATIENT_MAP` at it |
| `patient_alias_local.csv` | two columns, `clinical_name,patient_id_public` (for cases whose name is spelt differently in the clinical table and the cohort directory) | as above, or set `PAPER1_PATIENT_ALIAS` |

Both are listed in `.gitignore`. Without the first file the scripts fall back to using identifiers as they
are, which means they can only read released data that already uses the `P001`–`P109` codes.
