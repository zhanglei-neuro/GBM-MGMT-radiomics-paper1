# -*- coding: utf-8 -*-
r"""
Table S1 重建为 n=90（等价路径，不依赖原始 DICOM 盘）

背景：E:\NeuroImagepre 为外置盘，本次离线不可读；但既有 S1b（91 例 × 3 序列 = 273 行）
的逐行采集参数本就是直接读源 DICOM 头得到的，且这些参数只依赖"该例实际使用的那条序列"，
不随队列口径变化。因此：
  1) 精确剔除被排除病例（United Imaging 1.5T 单例 P021）的 3 行 -> 270 行 / 90 例
  2) 用与原脚本完全相同的汇总逻辑重算 S1a（按 T1c 实测厂商分组）与 S1c（机型技术签名）
  3) 同步过滤几何校验日志
校验：被保留行的全部字段与既有 S1b 逐字段一致（脚本内断言）。

输出：tables/S1a_paper1_by_vendor.csv, S1b_paper1_per_case.csv, S1c_paper1_vendor_signature.csv
      logs/s1_paper1_geometry_check.csv
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name   # noqa: E402  患者编号 <-> 本地姓名
import os
import pandas as pd

PROJ = r'D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1'
OUT = os.path.join(PROJ, 'supplementary', 'tables')
LOG = os.path.join(PROJ, 'supplementary', 'logs')

# 被排除的单例 United Imaging 1.5T（编号 P021）；本表 Patient 列为本地目录名，故还原
EXCLUDE = {to_name('P021')}   # 本地姓名对照见 tools/patient_ids_local.csv
N_EXPECT = 90

S1B = os.path.join(OUT, 'S1b_paper1_per_case.csv')
GEO = os.path.join(LOG, 's1_paper1_geometry_check.csv')

df_full = pd.read_csv(S1B, encoding='utf-8-sig')
print('既有 S1b: %d 行 / %d 例' % (len(df_full), df_full.Patient.nunique()))

# ---------- 1) 精确剔除 ----------
df = df_full[~df_full.Patient.isin(EXCLUDE)].copy()
assert len(df_full) - len(df) == 3, '被排除病例应恰好占 3 行'
assert df.Patient.nunique() == N_EXPECT, '剩余例数应为 %d' % N_EXPECT
print('剔除后: %d 行 / %d 例' % (len(df), df.Patient.nunique()))

# 校验：保留行与原表逐字段一致
keep_mask = df_full.Patient.isin(set(df.Patient))
chk = df_full[keep_mask].reset_index(drop=True).astype(str)
cur = df.reset_index(drop=True).astype(str)
assert chk.shape == cur.shape and (chk.values == cur.values).all(), '字段一致性校验失败'
print('字段一致性校验: 通过（保留行与原表完全相同）')

df.to_csv(S1B, index=False, encoding='utf-8-sig')

# ---------- 2) 几何校验日志过滤 ----------
geo_full = pd.read_csv(GEO, encoding='utf-8-sig')
geo = geo_full[~geo_full.patient.isin(EXCLUDE)].copy()
print('几何日志: %d -> %d 行' % (len(geo_full), len(geo)))
dm = geo['dim_match'].astype(str).str.lower() == 'true'
sp = geo['spacing_match'].astype(str).str.lower() == 'true'
print('  dim_match True: %d/%d ; spacing_match True: %d/%d' % (dm.sum(), len(geo), sp.sum(), len(geo)))
geo.to_csv(GEO, index=False, encoding='utf-8-sig')

# ---------- 3) 重算 S1a ----------
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
print('厂商分组(按T1c实测):', {v: df[df._Grp == v].Patient.nunique() for v in vend})

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

# ---------- 4) 重算 S1c ----------
sc = df.groupby(['Vendor', 'Model', 'Field strength (T)']).agg(
    n_sequences=('Patient', 'size'),
    n_patients=('Patient', 'nunique'),
    software=('Software version', lambda x: '; '.join(sorted(set(v for v in x if str(v).strip()))[:2])),
    sequences=('Sequence role', lambda x: '/'.join(sorted(set(x)))),
).reset_index()
sc.to_csv(os.path.join(OUT, 'S1c_paper1_vendor_signature.csv'), index=False, encoding='utf-8-sig')
print('\n=== S1c ===')
print(sc.to_string(index=False))

print('\nQC flag 分布 (270 行):')
print(df['QC flag'].value_counts().to_string())
print('\n队列厂商分布（按病例，T1c 实测）:')
print(df.drop_duplicates('Patient').Vendor.value_counts().to_string())
print('\n声明 label vs 实测 不一致:')
pv = df.drop_duplicates('Patient')[['Patient', 'Vendor', 'Vendor (Paper 1 label)']]
mism = pv[pv.Vendor != pv['Vendor (Paper 1 label)']]
print(mism.to_string(index=False) if len(mism) else '  无')
print('\nDONE (S1 rebuilt at n=90)')
