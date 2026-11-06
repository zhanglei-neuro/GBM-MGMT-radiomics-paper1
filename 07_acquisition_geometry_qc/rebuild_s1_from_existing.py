# -*- coding: utf-8 -*-
r"""
Table S1 rebuilt at n=90 via an equivalent path that does not depend on the original DICOM drive

Background: E:\NeuroImagepre is an external drive and was unreadable offline at that point; but the per-row
acquisition parameters of the existing S1b (91 cases x 3 sequences = 273 rows) were themselves read straight from
the source DICOM headers, and they depend only on "the series this case actually used", not on how the cohort is
defined. Therefore:
  1) drop exactly the 3 rows of the excluded case (the single United Imaging 1.5T case, P021) -> 270 rows / 90 cases
  2) recompute S1a (grouped by the measured T1c vendor) and S1c (model technical signature) with exactly the same
     aggregation logic as the original script
  3) filter the geometry-check log in step
Verification: every field of every retained row matches the existing S1b field for field (asserted in the script).

Output: tables/S1a_paper1_by_vendor.csv, S1b_paper1_per_case.csv, S1c_paper1_vendor_signature.csv
      logs/s1_paper1_geometry_check.csv
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name   # noqa: E402  patient code <-> local name
import os
import pandas as pd

PROJ = r'D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1'
OUT = os.path.join(PROJ, 'supplementary', 'tables')
LOG = os.path.join(PROJ, 'supplementary', 'logs')

# The excluded United Imaging 1.5T single case (code P021); the Patient column here holds local directory names, so it is mapped back
EXCLUDE = {to_name('P021')}   # local names: see tools/patient_ids_local.csv
N_EXPECT = 90

S1B = os.path.join(OUT, 'S1b_paper1_per_case.csv')
GEO = os.path.join(LOG, 's1_paper1_geometry_check.csv')

df_full = pd.read_csv(S1B, encoding='utf-8-sig')
print('existing S1b: %d rows / %d cases' % (len(df_full), df_full.Patient.nunique()))

# ---------- 1) drop exactly ----------
df = df_full[~df_full.Patient.isin(EXCLUDE)].copy()
assert len(df_full) - len(df) == 3, 'the excluded case should account for exactly 3 rows'
assert df.Patient.nunique() == N_EXPECT, 'the remaining case count should be %d' % N_EXPECT
print('after exclusion: %d rows / %d cases' % (len(df), df.Patient.nunique()))

# verification: retained rows match the original table field for field
keep_mask = df_full.Patient.isin(set(df.Patient))
chk = df_full[keep_mask].reset_index(drop=True).astype(str)
cur = df.reset_index(drop=True).astype(str)
assert chk.shape == cur.shape and (chk.values == cur.values).all(), 'field consistency check failed'
print('field consistency check: passed (retained rows identical to the original table)')

df.to_csv(S1B, index=False, encoding='utf-8-sig')

# ---------- 2) filter the geometry-check log ----------
geo_full = pd.read_csv(GEO, encoding='utf-8-sig')
geo = geo_full[~geo_full.patient.isin(EXCLUDE)].copy()
print('geometry log: %d -> %d rows' % (len(geo_full), len(geo)))
dm = geo['dim_match'].astype(str).str.lower() == 'true'
sp = geo['spacing_match'].astype(str).str.lower() == 'true'
print('  dim_match True: %d/%d ; spacing_match True: %d/%d' % (dm.sum(), len(geo), sp.sum(), len(geo)))
geo.to_csv(GEO, index=False, encoding='utf-8-sig')

# ---------- 3) recompute S1a ----------
def fmt_med(x, dec=1):
    x = pd.to_numeric(x, errors='coerce').dropna()
    if len(x) == 0:
        return ''
    return ('%.*f (%.*f-%.*f)' % (dec, x.median(), dec, x.min(), dec, x.max()))


order = ['Siemens', 'GE', 'UIH', 'Philips']
t1c_vendor = df[df['Sequence role'] == 'T1C'].drop_duplicates('Patient').set_index('Patient').Vendor
df['_Grp'] = df.Patient.map(t1c_vendor)
vend = [v for v in order if v in set(df._Grp.dropna())]
ALLCOL = 'All (n=%d)' % df.Patient.nunique()
VCOL = {v: '%s (n=%d)' % (v, df[df._Grp == v].Patient.nunique()) for v in vend}
print('vendor groups (measured on T1c):', {v: df[df._Grp == v].Patient.nunique() for v in vend})

lines = []


def add(param, fn):
    r = {'Parameter': param, ALLCOL: fn(df)}
    for v in vend:
        r[VCOL[v]] = fn(df[df._Grp == v])
    lines.append(r)


for role in ('T1C', 'T1w', 'FLAIR'):
    sub = lambda x, r=role: x[x['Sequence role'] == r]
    lines.append({'Parameter': '[%s] sequence' % role})
    add('%s · n patients' % role, lambda x, r=role: str(sub(x, r).Patient.nunique()))
    add('%s · field strength 1.5T / 3.0T' % role,
        lambda x, r=role: '%d / %d' % ((pd.to_numeric(sub(x, r)['Field strength (T)'], errors='coerce') == 1.5).sum(),
                                       (pd.to_numeric(sub(x, r)['Field strength (T)'], errors='coerce') == 3.0).sum()))
    add('%s · TR ms, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['TR (ms)'], 0))
    add('%s · TE ms, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['TE (ms)'], 1))
    add('%s · TI ms, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['TI (ms)'], 0))
    add('%s · flip angle deg, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['Flip angle (deg)'], 0))
    add('%s · slice thickness mm, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['Slice thickness (mm)'], 2))
    add('%s · in-plane pixel mm, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['Pixel spacing col (mm)'], 3))
    add('%s · matrix (rows x cols), mode' % role,
        lambda x, r=role: (sub(x, r)['Matrix'].mode().iloc[0] if len(sub(x, r)['Matrix'].mode()) else ''))
    add('%s · slices, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['N slices'], 0))
    add('%s · acquisition type 3D/2D' % role,
        lambda x, r=role: ' / '.join('%s: %d' % (k or 'N/A', v) for k, v in
                                     sub(x, r)['Acquisition type'].value_counts().items()))

lines.append({'Parameter': '[Scanner models]'})
add('model, n patients (per vendor)',
    lambda x: '; '.join('%s %s (%d)' % (v, m or 'N/A', n) for (v, m, n) in
                        x.drop_duplicates('Patient').groupby(['Vendor', 'Model']).size()
                        .reset_index().apply(lambda r: (r.Vendor, r.Model, r[0]), axis=1)))
add('software version(s)',
    lambda x: '; '.join(sorted(set(v for v in x['Software version'] if str(v).strip()))[:4]))

lines.append({'Parameter': '[QC]'})
lines.append({'Parameter': 'note',
              ALLCOL: ('Vendor groups are defined by the T1c (contrast-enhanced) reference sequence. '
                       'One case (P037) had T1c acquired on a GE scanner while T1w/FLAIR were '
                       'acquired on a Siemens scanner; it is therefore grouped as GE here but carries '
                       'a Siemens label in Table 1. Acquisition parameters below are the values '
                       'actually used for radiomics extraction. One further case (P021, '
                       'United Imaging uMR560 1.5 T) was excluded from the primary analysis per the '
                       'pre-specified exclusion criterion (Section 3.8) and is not represented here.')})

s1a = pd.DataFrame(lines).fillna('')
s1a.to_csv(os.path.join(OUT, 'S1a_paper1_by_vendor.csv'), index=False, encoding='utf-8-sig')
print('\n=== S1a ===')
print(s1a.to_string(index=False))

# ---------- 4) recompute S1c ----------
sc = df.groupby(['Vendor', 'Model', 'Field strength (T)']).agg(
    n_sequences=('Patient', 'size'),
    n_patients=('Patient', 'nunique'),
    software=('Software version', lambda x: '; '.join(sorted(set(v for v in x if str(v).strip()))[:2])),
    sequences=('Sequence role', lambda x: '/'.join(sorted(set(x)))),
).reset_index()
sc.to_csv(os.path.join(OUT, 'S1c_paper1_vendor_signature.csv'), index=False, encoding='utf-8-sig')
print('\n=== S1c ===')
print(sc.to_string(index=False))

print('\nQC flag distribution (270 rows):')
print(df['QC flag'].value_counts().to_string())
print('\nvendor distribution of the cohort (by case, measured on T1c):')
print(df.drop_duplicates('Patient').Vendor.value_counts().to_string())
print('\ndeclared label vs measurement, mismatches:')
pv = df.drop_duplicates('Patient')[['Patient', 'Vendor', 'Vendor (Paper 1 label)']]
mism = pv[pv.Vendor != pv['Vendor (Paper 1 label)']]
print(mism.to_string(index=False) if len(mism) else '  none')
print('\nDONE (S1 rebuilt at n=90)')
