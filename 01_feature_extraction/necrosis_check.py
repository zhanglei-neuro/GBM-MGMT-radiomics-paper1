# -*- coding: utf-8 -*-
r"""Check whether the Necrosis masks are empty (label values and the volume distribution of each subregion).

Usage:
    python necrosis_check.py <segmentation_tree> [report.txt]

<segmentation_tree> is the case level of a Raidionics output root, for example D:/Raidionics_data/patients/2
(one subdirectory per case, each holding T0/). Use it after a batch to decide whether the Necrosis channel
silently wrote empty masks.
"""
import os, sys, io

import numpy as np
import SimpleITK as sitk

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PAT = sys.argv[1] if len(sys.argv) > 1 else r'D:/Raidionics_data/patients/2'
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.getcwd(), 'necrosis_chk.txt')

pats = [d for d in sorted(os.listdir(PAT)) if os.path.isdir(os.path.join(PAT, d))]

lines = []
lines.append('%-12s %-12s %-22s %-22s' % ('patient', 'mask', 'unique values', 'non-zero voxels'))
lines.append('-' * 70)

rows = []
for p in pats:
    t0 = os.path.join(PAT, p, 'T0')
    for k, suf in (('TumorCore', 'annotation-Tumor_MRI_TumorCore.nii.gz'),
                   ('Necrosis', 'annotation-Necrosis_MRI_Necrosis.nii.gz'),
                   ('Edema', 'annotation-Edema.nii.gz'),
                   ('FLAIRChg', 'annotation-FLAIRChanges_MRI_FLAIRChanges.nii.gz')):
        mf = os.path.join(t0, '%s_t1-gd_pre_%s' % (p.lower(), suf))
        if not os.path.exists(mf):
            mf = os.path.join(t0, '%s_t2-flair_pre_%s' % (p.lower(), suf))
        if not os.path.exists(mf):
            rows.append((p, k, None, None)); continue
        a = sitk.GetArrayFromImage(sitk.ReadImage(mf))
        u = np.unique(a)
        nz = int((a > 0).sum())
        rows.append((p, k, u[:6].tolist(), nz))

# summary
for k in ('TumorCore', 'Necrosis', 'Edema', 'FLAIRChg'):
    sub = [r for r in rows if r[1] == k]
    empty = [r[0] for r in sub if r[3] == 0]
    ok = [r for r in sub if r[3]]
    vals = sorted(set(tuple(r[2]) for r in sub if r[2] is not None))
    lines.append('\n### %s : %d cases, %d empty, %d non-empty' % (k, len(sub), len(empty), len(ok)))
    lines.append('    label values seen: %s' % vals[:6])
    if empty:
        lines.append('    cases with an empty mask (first 30): %s' % empty[:30])

lines.append('\n\n### case-by-case detail (first 12)')
for r in rows[:48]:
    lines.append('%-12s %-12s %-22s %-22s' % (r[0], r[1], r[2], r[3]))

open(OUT, 'w', encoding='utf-8').write('\n'.join(lines))
print('written', OUT)
