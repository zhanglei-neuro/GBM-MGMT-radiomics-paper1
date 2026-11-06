# -*- coding: utf-8 -*-
r"""Compute the subregion volumes (mL) for every case -- used by Figure S2 and Table S3.
Cohort = the patients listed in analysis/predictions.csv (i.e. the Paper 1 primary analysis cohort, n=89).
Output: analysis/per_case_volumes.csv
"""
import os
import numpy as np
import pandas as pd
import nibabel as nib

import sys as _sys
_sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id, to_name   # noqa: E402  patient ID <-> local name

BASE = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
MAIN = r'D:/Raidionics_data/patients/2'          # TumorCore / Edema / FLAIRChanges / Brain
NEC = r'D:/Raidionics_data/patients_nec2/2'      # Necrosis (after repair)


def find(dd, kw):
    if os.path.isdir(dd):
        for f in sorted(os.listdir(dd)):
            if kw.lower() in f.lower():
                return os.path.join(dd, f)
    return None


pred = pd.read_csv(os.path.join(BASE, 'analysis', 'predictions.csv'), encoding='utf-8-sig')
pts = sorted(str(to_id(p)) for p in pred.patient.astype(str).tolist())   # output keys are always Pxxx
print('cohort from predictions.csv: n =', len(pts))

rows = []
for pt in pts:
    rec = {'patient': pt}
    dn = to_name(pt)          # data keys are Pxxx, image directories use the local name
    for k, dd, kw in [('TC', MAIN, 'annotation-Tumor_MRI_TumorCore'),
                      ('ED', MAIN, 'annotation-Edema'),
                      ('FL', MAIN, 'annotation-FLAIRChanges'),
                      ('BR', MAIN, 'annotation-Brain_MRI_Brain'),
                      ('NE', NEC, 'annotation-Necrosis')]:
        f = find(os.path.join(dd, dn, 'T0'), kw)
        if f is None:
            rec[k] = np.nan
            continue
        im = nib.load(f)
        z = im.header.get_zooms()
        rec[k] = float((np.asanyarray(im.dataobj) > 0).sum() * z[0] * z[1] * z[2] / 1000.0)
    rows.append(rec)

v = pd.DataFrame(rows)
out = os.path.join(BASE, 'analysis', 'per_case_volumes.csv')
v.to_csv(out, index=False, encoding='utf-8-sig')
print(v.describe().round(2).to_string())
print('missing:', v.isna().sum().to_dict())
print('saved ->', out)
