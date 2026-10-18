# -*- coding: utf-8 -*-
r"""
Supplementary S3 / Figure S2 -- 采集几何（层厚）对 Raidionics 分割的影响
A) 矢状位重建（同一物理 mm 视野）：3D 各向同性 (1 mm) vs 2D 厚层 (5 mm)，叠加亚区掩膜
B) 全脑掩膜体积按采集几何分组（OLS 校正年龄/性别）
C) 坏死灶 / (肿瘤核心+坏死) 比例（部分容积最敏感指标）
D) z 向覆盖 vs 层厚
输出: supplementary/figures/FigureS2_slice_thickness.png
      supplementary/tables/S3_slice_thickness.csv (+ .md)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name, to_id   # noqa: E402  患者编号 <-> 本地姓名
import os
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
import statsmodels.api as sm

BASE = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
SUP = os.path.join(BASE, 'supplementary')
FIGD = os.path.join(SUP, 'figures'); os.makedirs(FIGD, exist_ok=True)
INP = r'D:/Raidionics_data/input/2'
MAIN = r'D:/Raidionics_data/patients/2'
NEC = r'D:/Raidionics_data/patients_nec2/2'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 120})

# 键为患者目录名（以姓名命名），由去标识编号经 to_name() 还原
CASES = [(to_name('P098'), 'Siemens Verio 3.0 T | 3D MP-RAGE\n1.0 x 1.0 x 1.0 mm', 'thin'),
         (to_name('P010'), 'Siemens Verio 3.0 T | 2D T1-TIRM\n0.34 x 0.34 x 5.0 mm (6 mm spacing)', 'thick')]


def find(d, kw):
    if not os.path.isdir(d):
        return None
    for f in sorted(os.listdir(d)):
        if kw.lower() in f.lower():
            return os.path.join(d, f)
    return None


def get(pat):
    img = nib.as_closest_canonical(nib.load(find(os.path.join(INP, pat, 'T0'), 't1-gd')))
    out = {'img': img}
    for k, dd, kw in [('TC', MAIN, 'annotation-Tumor_MRI_TumorCore'),
                      ('ED', MAIN, 'annotation-Edema'),
                      ('NE', NEC, 'annotation-Necrosis')]:
        f = find(os.path.join(dd, pat, 'T0'), kw)
        out[k] = nib.as_closest_canonical(nib.load(f)) if f else None
    return out


fig = plt.figure(figsize=(10.2, 6.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], hspace=0.40, wspace=0.32)
YW, ZW = 95.0, 62.0      # 共同物理视野 (mm): y ±95, z ±62

for i, (pat, lab, grp) in enumerate(CASES):
    D = get(pat)
    img = D['img']
    vol = np.asanyarray(img.dataobj).astype(float)
    px, py, pz = [float(x) for x in img.header.get_zooms()[:3]]
    tc = np.asanyarray(D['TC'].dataobj) > 0
    ne = np.asanyarray(D['NE'].dataobj) > 0 if D['NE'] is not None else np.zeros_like(tc)
    ed = np.asanyarray(D['ED'].dataobj) > 0 if D['ED'] is not None else np.zeros_like(tc)
    ci = int(np.mean(np.where(tc.any(axis=(1, 2)))[0]))   # R-L index
    cj = int(np.mean(np.where(tc.any(axis=(0, 2)))[0]))   # A-P index
    ck = int(np.mean(np.where(tc.any(axis=(0, 1)))[0]))   # S-I index
    print(pat, 'shape', tc.shape, 'centroid', (ci, cj, ck), 'zooms', (px, py, pz))

    def crop(vol3, cj=cj, ck=ck, py=py, pz=pz, ci=ci):
        sl = vol3[ci]
        ym = (np.arange(sl.shape[0]) - cj) * py
        zm = (np.arange(sl.shape[1]) - ck) * pz
        ys = np.abs(ym) <= YW; zs = np.abs(zm) <= ZW
        return sl[np.ix_(ys, zs)].T[::-1], ym[ys], zm[zs]

    sub, ym, zm = crop(vol)
    ext = [ym.min(), ym.max(), zm.max(), zm.min()]
    lo, hi = np.percentile(sub, 2), np.percentile(sub, 99.5)
    ax = fig.add_subplot(gs[0, i])
    ax.imshow(np.clip(sub, lo, hi), cmap='gray', extent=ext, origin='upper',
              aspect='equal', interpolation='nearest')
    ov = np.zeros((*sub.shape, 4))
    for msk, col in [(ed, (0.20, 0.55, 0.95)), (tc, (0.90, 0.20, 0.20)), (ne, (0.15, 0.75, 0.30))]:
        mm, _, _ = crop(msk.astype(float))
        ov[mm > 0] = [*col, 0.50]
    ax.imshow(ov, extent=ext, origin='upper', aspect='equal', interpolation='nearest')
    ax.set_xlim(-YW, YW); ax.set_ylim(-ZW, ZW)
    ax.set_xlabel('A-P (mm)', fontsize=8); ax.set_ylabel('S-I (mm)', fontsize=8) if i == 0 else None
    ax.tick_params(labelsize=7)
    ax.set_title('%s\n%s' % (pat, lab), fontsize=8.2)
    ax.text(0.02, 0.02, 'voxel %.2f x %.2f x %.2f mm' % (px, py, pz), transform=ax.transAxes,
            fontsize=7, bbox=dict(fc='w', alpha=0.8, ec='none', pad=1.5))
    if grp == 'thick':
        for zz in np.arange(zm.min(), zm.max(), pz):
            ax.axhline(zz, color='yellow', lw=0.5, alpha=0.55)

axL = fig.add_subplot(gs[0, 2]); axL.axis('off')
axL.legend(handles=[Patch(fc=(0.20, 0.55, 0.95), alpha=.6, label='Edema / SNFH'),
                    Patch(fc=(0.90, 0.20, 0.20), alpha=.6, label='Tumour core'),
                    Patch(fc=(0.15, 0.75, 0.30), alpha=.6, label='Necrosis')],
           loc='upper center', frameon=False, fontsize=9)
axL.text(0.5, 0.58,
         'Both panels: sagittal reformats at the\n'
         'same physical scale (mm).\n\n'
         'Yellow lines = native slice planes of\nthe 2D acquisition (6 mm apart).\n\n'
         'Raidionics zero-fills every volume to\n1 mm isotropic before inference\n'
         '(output_spacing = 1,1,1), as Zahirovic\net al. did via 1 mm template registration.\n\n'
         'But interpolation cannot restore detail\n'
         'that was never sampled: effective\n'
         'through-plane resolution stays 5 mm,\n'
         'and 1 mm interslice gaps (17% of\nz-coverage) are never acquired.',
         ha='center', va='top', fontsize=7.2, linespacing=1.45)

# ---------------- B / C / D ----------------
v = pd.read_csv(os.path.join(BASE, 'analysis/per_case_volumes.csv'))
v['patient'] = [to_id(p) for p in v['patient']]     # 落盘表统一归一化为 Pxxx
s1b = pd.read_csv(os.path.join(SUP, 'tables/S1b_paper1_per_case.csv'))
s1b['Patient'] = [to_id(p) for p in s1b['Patient']]   # 落盘表统一归一化为 Pxxx
t = s1b[s1b['Sequence role'] == 'T1C'].copy()
t['th'] = pd.to_numeric(t['Slice thickness (mm)'], errors='coerce')
geo = pd.read_csv(os.path.join(SUP, 'logs/s1_paper1_geometry_check.csv'))
geo = geo[geo.role == 'T1C'].copy()
geo['dim'] = geo.nifti_dim.astype(str)
geo['nz'] = geo.dim.str.split('x').apply(lambda s: min(int(x) for x in s))
geo['zmax_mm'] = geo.nifti_zoom.astype(str).str.split('x').apply(lambda s: max(float(x) for x in s))
geo['cov'] = geo.nz * geo.zmax_mm
m = v.merge(t[['Patient', 'th', 'Vendor', 'Acquisition type', 'N slices']],
            left_on='patient', right_on='Patient', how='left')
m = m.merge(geo[['patient', 'cov']], on='patient', how='left')
m['grp'] = np.where(m.th >= 3, '2D thick\n(5-6 mm)', '3D thin\n(1.0-1.2 mm)')
m['NE_frac'] = m.NE / (m.TC + m.NE).replace(0, np.nan)
co = pd.read_csv(r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts/cohort_MGMT.csv',
                 encoding='utf-8-sig').set_index('patient')
m['age'] = pd.to_numeric(co.loc[m.patient, 'age'].values, errors='coerce')
m['male'] = co.loc[m.patient, 'gender'].astype(str).str.strip().eq('M').astype(int).values
G = ['3D thin\n(1.0-1.2 mm)', '2D thick\n(5-6 mm)']
COL = ['#3b6ea5', '#c0504d']

axB = fig.add_subplot(gs[1, 0])
data = [m.loc[m.grp == g, 'BR'].dropna() / 1000.0 for g in G]
bp = axB.boxplot(data, widths=.5, patch_artist=True, showfliers=False,
                 medianprops=dict(color='k', lw=1.4))
for b, c in zip(bp['boxes'], COL):
    b.set_facecolor(c); b.set_alpha(.35)
for i, d in enumerate(data):
    axB.scatter(np.random.default_rng(1).normal(i + 1, .07, len(d)), d, s=9, color=COL[i], alpha=.75, lw=0)
ok = m[['BR', 'th', 'age', 'male']].dropna()
fit = sm.OLS(ok.BR.values / 1000.0,
             sm.add_constant(pd.DataFrame({'thick': (ok.th >= 3).astype(int),
                                           'age': ok.age, 'male': ok.male}))).fit()
axB.set_xticklabels(G, fontsize=8); axB.set_ylabel('Brain mask volume (L)', fontsize=8.5)
axB.set_title('B. Brain coverage truncated\n(adj. $\\Delta$=%.2f L, p=%.3f)' %
              (fit.params['thick'], fit.pvalues['thick']), fontsize=8.5)

axC = fig.add_subplot(gs[1, 1])
data = [m.loc[m.grp == g, 'NE_frac'].dropna() for g in G]
bp = axC.boxplot(data, widths=.5, patch_artist=True, showfliers=False,
                 medianprops=dict(color='k', lw=1.4))
for b, c in zip(bp['boxes'], COL):
    b.set_facecolor(c); b.set_alpha(.35)
for i, d in enumerate(data):
    axC.scatter(np.random.default_rng(2).normal(i + 1, .07, len(d)), d, s=9, color=COL[i], alpha=.75, lw=0)
axC.set_xticklabels(G, fontsize=8); axC.set_ylabel('Necrosis / (TC + necrosis)', fontsize=8.5)
axC.set_title('C. No bias in the most PV-sensitive\nratio (p=%.2f)' %
              stats.mannwhitneyu(data[0], data[1]).pvalue, fontsize=8.5)

axD = fig.add_subplot(gs[1, 2])
for g, c, nm in [(G[0], COL[0], '3D thin'), (G[1], COL[1], '2D thick')]:
    d = m[m.grp == g]
    axD.scatter(d['th'], d['cov'], s=22, color=c, alpha=.85, lw=0, label='%s (n=%d)' % (nm, len(d)))
axD.axhspan(140, 150, color='gray', alpha=.25, lw=0)
axD.text(1.6, 152, 'typical brain S-I extent (140-150 mm)', fontsize=6.8, color='dimgray')
axD.set_xlabel('T1c slice thickness (mm)', fontsize=8.5)
axD.set_ylabel('z-coverage (mm)', fontsize=8.5)
axD.legend(frameon=False, fontsize=7.2, loc='center right')
axD.set_title('D. z-coverage vs slice thickness', fontsize=8.5)

fig.suptitle('Figure S2. Acquisition geometry (slice thickness) and its effect on automatic segmentation',
             fontsize=10.5, y=0.985)
out = os.path.join(FIGD, 'FigureS2_slice_thickness.png')
fig.savefig(out, bbox_inches='tight'); plt.close(fig)
print('saved', out)

# ---------------- 表 S3 ----------------
def iqr(x):
    return '%.2f [%.2f-%.2f]' % (x.median(), x.quantile(.25), x.quantile(.75))


rows = []
for c, nm in [('BR', 'Whole-brain mask volume (mL)'), ('TC', 'Tumour core volume (mL)'),
              ('NE', 'Necrosis volume (mL)'), ('ED', 'Edema volume (mL)'),
              ('FL', 'FLAIRChanges volume (mL)'), ('NE_frac', 'Necrosis / (TC+Necrosis)')]:
    a = m.loc[m.grp == G[0], c].dropna(); b = m.loc[m.grp == G[1], c].dropna()
    u = stats.mannwhitneyu(a, b)
    rows.append({'Metric': nm,
                 '3D thin (n=%d), median [IQR]' % len(a): iqr(a),
                 '2D thick (n=%d), median [IQR]' % len(b): iqr(b),
                 'p (Mann-Whitney)': round(u.pvalue, 4),
                 'Rank-biserial r': round(2 * u.statistic / (len(a) * len(b)) - 1, 3)})
t3 = pd.DataFrame(rows)
t3.to_csv(os.path.join(SUP, 'tables/S3_slice_thickness.csv'), index=False, encoding='utf-8-sig')
with open(os.path.join(SUP, 'tables/S3_slice_thickness.md'), 'w', encoding='utf-8') as f:
    f.write('### Table S3. Segmentation-derived metrics stratified by T1c acquisition geometry\n\n')
    f.write(t3.to_markdown(index=False))
    f.write('\n\n**Whole-brain mask volume -- OLS adjusted for age and sex (n=%d):**\n\n' % len(ok))
    f.write('| term | coef (L) | 95% CI | p |\n|---|---|---|---|\n')
    for tm_ in ['thick', 'age', 'male']:
        f.write('| %s | %.4f | %.4f to %.4f | %.4f |\n' % (
            tm_, fit.params[tm_], fit.conf_int().loc[tm_, 0], fit.conf_int().loc[tm_, 1], fit.pvalues[tm_]))
    f.write('\n*2D thick = slice thickness >= 5 mm (all 2D acquisitions); '
            '3D thin = <= 1.2 mm (all 3D acquisitions).*\n')
print(t3.to_string(index=False))
print('\nAdjusted effect of 2D thick on brain volume: %.4f L (95%%CI %.4f to %.4f), p=%.4f'
      % (fit.params['thick'], fit.conf_int().loc['thick', 0], fit.conf_int().loc['thick', 1],
         fit.pvalues['thick']))
