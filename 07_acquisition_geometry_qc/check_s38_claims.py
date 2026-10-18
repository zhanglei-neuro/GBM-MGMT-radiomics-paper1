# -*- coding: utf-8 -*-
r"""核验 3.8 节两个细节声明（n=89 主分析队列，3D 薄层 59 vs 2D 厚层 30）：
1) 小坏死灶（<0.1 mL）检出率两组是否相同
2) 脑掩膜触及扫描 z 向首层/末层边界的例数（按组）
输出: supplementary/scripts/s38_claims_check.txt
"""
import os
import glob
import numpy as np
import pandas as pd

import sys as _sys
_sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import SimpleITK as sitk

OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
MAIN = r'D:\Raidionics_data\patients\2'
NEC = r'D:\Raidionics_data\patients_nec2\2'
lines = []

s1b = pd.read_csv(os.path.join(OUT, 'supplementary/tables/S1b_paper1_per_case.csv'))
s1b['Patient'] = [to_id(p) for p in s1b['Patient']]
t1c = s1b[s1b['Sequence role'] == 'T1C'].copy()
t1c['th'] = pd.to_numeric(t1c['Slice thickness (mm)'], errors='coerce')
thin = set(t1c.loc[t1c.th <= 1.2, 'Patient'])
thick = set(t1c.loc[t1c.th >= 3, 'Patient'])

vol = pd.read_csv(os.path.join(OUT, 'analysis/per_case_volumes.csv'))
vol['patient'] = [to_id(p) for p in vol['patient']]
vol['grp'] = np.where(vol.patient.isin(thick), '2D_thick',
                      np.where(vol.patient.isin(thin), '3D_thin', 'other'))
vol = vol[vol.grp != 'other'].copy()
lines.append('n: %s' % vol.grp.value_counts().to_dict())

# ---- 声明 1: 小坏死灶 (<0.1 mL) 检出率 ----
# 1 mL = 1000 mm^3; per_case_volumes 的 NE 单位为 mL（与表 S3 一致）
vol['nec_small'] = vol['NE'] < 0.1
rate = vol.groupby('grp')['nec_small'].agg(['sum', 'count', 'mean'])
lines.append('\n[claim 1] small-necrosis (<0.1 mL) detection rate')
lines.append(rate.to_string())
for g in ['3D_thin', '2D_thick']:
    s = vol[vol.grp == g]
    lines.append('  %s: %d/%d = %.0f%%' % (g, s.nec_small.sum(), len(s), 100 * s.nec_small.mean()))

# ---- 声明 2: 脑掩膜触及 z 向首层/末层 ----
lines.append('\n[claim 2] brain mask touching first/last z slice')


def brain_touches_z(p):
    pat = os.path.join(MAIN, p, 'T0', '*annotation-Brain_MRI_Brain.nii.gz')
    fs = glob.glob(pat)
    if not fs:
        return 'no_mask'
    a = sitk.GetArrayFromImage(sitk.ReadImage(fs[0])) > 0   # (z, y, x)
    if a.sum() == 0:
        return 'empty'
    z0 = a[0, :, :].any()
    z1 = a[-1, :, :].any()
    return 'touch' if (z0 or z1) else 'ok'


vol['brain_z'] = vol.patient.map(brain_touches_z)
lines.append(vol.groupby('grp')['brain_z'].value_counts().unstack(fill_value=0).to_string())
for g in ['3D_thin', '2D_thick']:
    s = vol[vol.grp == g]
    lines.append('  %s: touch %d / %d' % (g, int((s.brain_z == 'touch').sum()), len(s)))

txt = '\n'.join(lines)
print(txt)
open(os.path.join(OUT, 'supplementary/scripts/s38_claims_check.txt'), 'w',
     encoding='utf-8').write(txt + '\n')
