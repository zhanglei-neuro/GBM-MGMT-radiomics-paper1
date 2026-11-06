# -*- coding: utf-8 -*-
r"""
patient_ids.py -- patient identifier normalisation

The public data (Zenodo deposit and journal supplementary tables) use the de-identified codes **P001-P109**
(see `zenodo/README_DATASET.md` and the supplementary `TABLE_LEGENDS.md`).
Data tables in the local working environment use patient names in pinyin. Put the code <-> local name mapping at

    github/tools/patient_ids_local.csv      columns: patient_id_public,patient_name_original

That file contains PHI and **must not be committed** (it is listed in `github/.gitignore`).
If the local data already uses codes (for example because you are working straight from the Zenodo release),
the file is not needed.

Usage
----
    from patient_ids import to_id, to_name, normalize, drop_uih

    to_id('P001')        -> 'P001'      already a code; idempotent
    to_id('<PATIENT>')   -> 'P001'      local name -> code (case-insensitive)
    to_name('P001')      -> '<PATIENT>' code -> local name (for building local paths)

Three conventions

1. **Data keys** (cohort membership, the keys of per-case mappings, case selection for figures) are always Pxxx:
   after loading any patient-indexed table call `normalize()` before doing anything else.
   That way the same code works on both local data (names) and released data (codes).
2. **File paths** must use the local directory names: DICOM and NIfTI directories are named after the patient,
   so build paths with `to_name(pid)`.
3. **Cohort exclusion is always decided by vendor** (see `drop_uih` / `UIH_MASK`), never by name:
   in the released data names have been replaced by codes, so excluding by name **fails silently** and changes n.
"""
import os
import csv

_HERE = os.path.dirname(os.path.abspath(__file__))

UIH = 'UIH'          # vendor label for United Imaging (values of the vendor column in cohort_MGMT.csv)

_LOCAL_CANDIDATES = [
    os.environ.get('PAPER1_PATIENT_MAP', ''),
    os.path.join(_HERE, 'patient_ids_local.csv'),
    os.path.join(_HERE, '..', '..', '_internal', 'patient_ids_local.csv'),
]

_ALIAS_CANDIDATES = [
    os.environ.get('PAPER1_PATIENT_ALIAS', ''),
    os.path.join(_HERE, 'patient_alias_local.csv'),
    os.path.join(_HERE, '..', '..', '_internal', 'patient_alias_local.csv'),
]


def load_map():
    """Return (id2name, name2id); both are empty dictionaries when no mapping file is found."""
    id2name, name2id = {}, {}
    for p in _LOCAL_CANDIDATES:
        if p and os.path.isfile(p):
            with open(p, encoding='utf-8-sig', newline='') as f:
                r = csv.reader(f)
                next(r, None)
                for row in r:
                    if len(row) >= 2 and row[0].strip() and row[1].strip():
                        pid, nm = row[0].strip(), row[1].strip()
                        id2name[pid.upper()] = nm
                        name2id[nm.lower()] = pid
            break
    return id2name, name2id


ID2NAME, NAME2ID = load_map()
HAS_MAP = bool(ID2NAME)


def to_id(x):
    """Name -> Pxxx; codes and unknown identifiers are returned unchanged."""
    s = str(x).strip()
    if not s:
        return x
    if s.upper() in ID2NAME:
        return s.upper()
    return NAME2ID.get(s.lower(), x)


def to_name(x):
    """Pxxx -> local name; anything that is not a code, or is unknown, is returned unchanged."""
    s = str(x).strip()
    return ID2NAME.get(s.upper(), x)


def load_alias():
    r"""
    Spelling variants in the clinical table -> Pxxx.

    A few cases are spelt differently in the Name column of `GBM_clinical_data.xlsx` than in the imaging cohort
    directory name (a different first syllable, or one letter short, for example). Put those pairs in the optional
    file

        github/tools/patient_alias_local.csv      columns: clinical_name,patient_id_public

    That file **contains no PHI** (only spelling variants and existing codes) and may be published with the
    repository; it is empty by default.
    """
    alias = {}
    for p in _ALIAS_CANDIDATES:
        if p and os.path.isfile(p):
            with open(p, encoding='utf-8-sig', newline='') as f:
                r = csv.reader(f)
                next(r, None)
                for row in r:
                    if len(row) >= 2 and row[0].strip() and row[1].strip():
                        alias[row[0].strip()] = row[1].strip()
                        alias[row[0].strip().lower()] = row[1].strip()
            break
    return alias


def resolve(x):
    """Look in the spelling-variant table first, then fall back to to_id; returns Pxxx (unknown identifiers unchanged)."""
    s = str(x).strip()
    if not s:
        return x
    a = load_alias()
    return to_id(a.get(s, a.get(s.lower(), s)))


def normalize(df, col=None, inplace=True):
    """Normalise the patient identifiers in a DataFrame index (or in column `col`) to Pxxx."""
    out = df if inplace else df.copy()
    if col is None:
        out.index = [to_id(p) for p in out.index]
    else:
        out[col] = [to_id(p) for p in out[col]]
    return out


def norm_seq(seq):
    """Normalise any sequence of patient identifiers to Pxxx."""
    return [to_id(x) for x in seq]


def drop_uih(df, col='vendor'):
    """Drop the single United Imaging case (decided from the vendor column, not by name)."""
    if df is None or col not in getattr(df, 'columns', []):
        return df
    return df[df[col].astype(str).str.strip() != UIH].copy()


def uih_mask(ids, co, col='vendor'):
    """Given a sequence of patient identifiers and a cohort table indexed by patient, return a boolean list marking UIH."""
    return [str(co.loc[p, col]).strip() == UIH for p in ids]
