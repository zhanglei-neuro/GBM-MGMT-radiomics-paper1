# -*- coding: utf-8 -*-
"""核查 Necrosis 掩码是否为空（标签取值 / 各亚区体积分布）。

用法:
    python necrosis_check.py <segmentation_tree> [report.txt]

<segmentation_tree> 为 Raidionics 输出根目录下的病例层，例如 D:/Raidionics_data/patients/2
（每例一个子目录，内含 T0/）。用于在主批跑后判定 Necrosis 通道是否静默输出了空掩码。
"""
import os, sys, io

import numpy as np
import SimpleITK as sitk

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

PAT = sys.argv[1] if len(sys.argv) > 1 else r'D:/Raidionics_data/patients/2'
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.getcwd(), 'necrosis_chk.txt')

pats = [d for d in sorted(os.listdir(PAT)) if os.path.isdir(os.path.join(PAT, d))]

lines = []
lines.append('%-12s %-12s %-22s %-22s' % ('患者', '掩码', '唯一值', '非零体素数'))
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

# 汇总
for k in ('TumorCore', 'Necrosis', 'Edema', 'FLAIRChg'):
    sub = [r for r in rows if r[1] == k]
    empty = [r[0] for r in sub if r[3] == 0]
    ok = [r for r in sub if r[3]]
    vals = sorted(set(tuple(r[2]) for r in sub if r[2] is not None))
    lines.append('\n### %s : 共 %d 例, 空掩码 %d 例, 非空 %d 例' % (k, len(sub), len(empty), len(ok)))
    lines.append('    出现过的标签取值: %s' % vals[:6])
    if empty:
        lines.append('    空掩码患者(前30): %s' % empty[:30])

lines.append('\n\n### 逐例明细 (前 12 例)')
for r in rows[:48]:
    lines.append('%-12s %-12s %-22s %-22s' % (r[0], r[1], r[2], r[3]))

open(OUT, 'w', encoding='utf-8').write('\n'.join(lines))
print('written', OUT)
