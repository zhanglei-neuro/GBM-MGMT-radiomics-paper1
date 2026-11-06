# -*- coding: utf-8 -*-
r"""Paper 1 - volumetric biomarker analysis (benchmarked against Zahirovic 2025, Front Neurosci)
v3 (2026-10-11): fixes the AUC label-polarity bug raised as M1 in peer review and switches to a more reliable volume source
  - v1 negated every ratio (score=-s), on the assumption, written in its own comments, that "MGMT+ has the lower ratio",
    which is the opposite of what was measured; that flipped the TT/ED and TC/ED AUCs to 1-AUC (0.383/0.367). From v2 on
    the natural direction is used throughout (score=+s, higher value -> predicts MGMT+), and the table states explicitly
    which group has the higher median.
  - v2 used PyRadiomics shape features as the volume source, which carried two spurious missing values (one TumorCore|T1CE
    extraction failure and one Edema mask that could not be transferred to the FLAIR grid); treating them as zero would
    introduce bias.
  - v3 uses analysis/per_case_volumes.csv: each subregion volume is counted directly on its own native grid (n=89, no
    missing values), and TumorCore/Necrosis/Edema share the T1c reference grid, so the ratios are geometrically consistent.
  - Added: within-family BH-FDR q values for the seven volumetric markers (peer-review Major-2 / M6) and
    Mann-Whitney p values for the per-vendor subgroups.
Outputs: analysis/metrics_volume.json, table3_volume.csv/md, figures/fig8_volume.png
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  patient code <-> local name
import os, json, warnings
import numpy as np, pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis'); FIG = os.path.join(OUT, 'figures')
os.makedirs(FIG, exist_ok=True)
SEED = 42
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 120,
                     'savefig.bbox': 'tight'})
def log(*a): print(*a, flush=True)

co = pd.read_csv(os.path.join(SRC, 'scripts', 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]            # all data keys are Pxxx
V = pd.read_csv(os.path.join(ANA, 'per_case_volumes.csv'), encoding='utf-8-sig').set_index('patient')
V.index = [to_id(p) for p in V.index]
common = [p for p in co.index if p in V.index]
# Exclude the single United Imaging 1.5T case (see manuscript 3.8). Decided by vendor rather than by name
# (in the released data patients are Pxxx, so excluding by name fails silently and changes n)
common = [p for p in common if str(co.loc[p, 'vendor']).strip() != 'UIH']
V = V.loc[common]; co = co.loc[common]
log('n=%d  volumes missing: %s' % (len(V), V.isna().sum().to_dict()))
y = co['y'].values; vendor = co['vendor'].values
age = pd.to_numeric(co['age'], errors='coerce').values

V = V.rename(columns={'TC': 'TumorCore', 'NE': 'Necrosis', 'ED': 'Edema',
                      'FL': 'FLAIRChanges', 'BR': 'Brain'})
V['TotalTumor'] = V['TumorCore'] + V['Necrosis']      # an empty necrosis mask means 0, which is biologically valid
V['TT_ED_ratio'] = V['TotalTumor'] / V['Edema'].replace(0, np.nan)
V['TC_ED_ratio'] = V['TumorCore'] / V['Edema'].replace(0, np.nan)

def auc_ci(score, yv, n=2000, seed=SEED):
    score = np.asarray(score, dtype=float); yv = np.asarray(yv)
    m = np.isfinite(score)
    score, yv = score[m], yv[m]
    a = roc_auc_score(yv, score)
    rng = np.random.default_rng(seed); bs = []
    for _ in range(n):
        i = rng.integers(0, len(yv), len(yv))
        if len(np.unique(yv[i])) < 2: continue
        bs.append(roc_auc_score(yv[i], score[i]))
    return a, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

MARKERS = ['TumorCore', 'Necrosis', 'Edema', 'FLAIRChanges', 'TotalTumor', 'TT_ED_ratio', 'TC_ED_ratio']
rows = []
log('\n%-16s %5s %-22s %-22s %8s %8s' % ('marker', 'n', 'MGMT- median[IQR]', 'MGMT+ median[IQR]', 'p', 'AUC'))
log('  (AUC in the natural direction: higher value -> predicts MGMT+; nothing negated or flipped)')
ps_raw = {}
for col in MARKERS:
    s = V[col].values; m = np.isfinite(s)
    if m.sum() < 20: continue
    a = s[m & (y == 0)]; b = s[m & (y == 1)]
    p = stats.mannwhitneyu(a, b).pvalue
    ps_raw[col] = p
    au, lo, hi = auc_ci(s[m], y[m])          # natural direction, no flipping
    higher_in = 'MGMT+' if np.median(b) > np.median(a) else 'MGMT-'
    q = lambda z: '%6.2f [%.2f-%.2f]' % (np.median(z), np.percentile(z, 25), np.percentile(z, 75))
    log('%-16s %5d %-22s %-22s %8.4f %8.3f' % (col, int(m.sum()), q(a), q(b), p, au))
    rows.append([col, int(m.sum()), higher_in, q(a), q(b), '%.4f' % p, '%.3f' % au, '%.3f–%.3f' % (lo, hi)])

# within-family BH-FDR (the seven volumetric markers form one prespecified family)
pser = pd.Series(ps_raw).sort_values()
k = len(pser)
qv = pser.values * k / (np.arange(k) + 1)
qv = np.minimum.accumulate(qv[::-1])[::-1]
qmap = dict(zip(pser.index, qv))
for r in rows:
    r.append('%.4f' % qmap[r[0]])
log('\nBH-FDR within the 7-marker family:')
for kk, vv2 in qmap.items():
    log('  %-14s p=%.4f  q=%.4f' % (kk, ps_raw[kk], vv2))

t3 = pd.DataFrame(rows, columns=['Marker', 'n', 'Higher in', 'MGMT-unmethylated median [IQR]',
                                 'MGMT-methylated median [IQR]', 'p (Mann–Whitney)',
                                 'AUC (predict MGMT+)', 'AUC 95% CI', 'BH-FDR q (7-marker family)'])
t3.to_csv(os.path.join(ANA, 'table3_volume.csv'), index=False, encoding='utf-8-sig')
open(os.path.join(ANA, 'table3_volume.md'), 'w', encoding='utf-8').write(t3.to_markdown(index=False))

# multivariable logistic: log(ratio) + age  (nested CV)
log('\n=== multivariate: log(TT/ED ratio) + age ===')
ratio = V['TT_ED_ratio'].values
Xm = pd.DataFrame({'log_ratio': np.log(ratio), 'age': age}, index=V.index)
mv = {}
skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
p = np.zeros(len(y)); coef_lr = []
for tr, te in skf.split(Xm, y):
    pipe = Pipeline([('i', SimpleImputer(strategy='median')), ('s', StandardScaler()),
                     ('m', LogisticRegression(max_iter=2000))])
    pipe.fit(Xm.iloc[tr], y[tr])
    p[te] = pipe.predict_proba(Xm.iloc[te])[:, 1]
    coef_lr.append(pipe.named_steps['m'].coef_[0])
mv['log_ratio_age_auc'] = float(roc_auc_score(y, p))
mv['mean_coef_log_ratio'] = float(np.mean([c[0] for c in coef_lr]))
mv['mean_coef_age'] = float(np.mean([c[1] for c in coef_lr]))
log('  AUC=%.3f  coef(log_ratio)=%.3f  coef(age)=%.3f' %
    (mv['log_ratio_age_auc'], mv['mean_coef_log_ratio'], mv['mean_coef_age']))

# univariate ratio only (natural direction)
au_r, lo_r, hi_r = auc_ci(V['TT_ED_ratio'].values, y)
mv['tt_ed_ratio_auc'] = au_r; mv['tt_ed_ratio_ci'] = [lo_r, hi_r]
mv['tt_ed_ratio_n'] = int(np.isfinite(ratio).sum())
mv['bh_q'] = {k2: float(v2) for k2, v2 in qmap.items()}
mv['p_raw'] = {k2: float(v2) for k2, v2 in ps_raw.items()}
mv['volume_source'] = 'analysis/per_case_volumes.csv (native-grid voxel counting, n=89)'

# per-vendor subgroups
log('\n=== univariate TT/ED ratio AUC by vendor (natural direction) ===')
sub = {}
for v in ['Siemens', 'GE']:
    m = (vendor == v) & np.isfinite(ratio)
    if m.sum() >= 15 and len(np.unique(y[m])) > 1:
        a, lo, hi = auc_ci(ratio[m], y[m])
        pv = stats.mannwhitneyu(ratio[m & (y == 0)], ratio[m & (y == 1)]).pvalue
        sub[v] = {'n': int(m.sum()), 'n_pos': int(y[m].sum()), 'auc': a, 'ci': [lo, hi], 'p': float(pv)}
        log('  %-8s n=%2d (+%d) AUC=%.3f (%.3f-%.3f)  p=%.4f' % (v, m.sum(), y[m].sum(), a, lo, hi, pv))
mv['by_vendor_tt_ed'] = sub
# TumorCore volume AUC by vendor (supports the "consistent direction" statement in the Discussion)
tc_sub = {}
for v in ['Siemens', 'GE']:
    m = (vendor == v) & np.isfinite(V['TumorCore'].values)
    if m.sum() >= 15 and len(np.unique(y[m])) > 1:
        a, lo, hi = auc_ci(V['TumorCore'].values[m], y[m])
        pv = stats.mannwhitneyu(V['TumorCore'].values[m & (y == 0)], V['TumorCore'].values[m & (y == 1)]).pvalue
        tc_sub[v] = {'n': int(m.sum()), 'n_pos': int(y[m].sum()), 'auc': a, 'ci': [lo, hi], 'p': float(pv)}
        log('  TumorCore %-8s n=%2d AUC=%.3f (%.3f-%.3f) p=%.4f' % (v, m.sum(), a, lo, hi, pv))
mv['by_vendor_tumorcore'] = tc_sub
json.dump(mv, open(os.path.join(ANA, 'metrics_volume.json'), 'w'), indent=2)

# ---------------- fig8 ----------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
# (a) ratio distributions
for lab, c, nm in [(0, '#999999', 'MGMT-unmethylated'), (1, '#C44E52', 'MGMT-methylated')]:
    d = ratio[y == lab]; d = d[np.isfinite(d)]
    ax[0].hist(d, bins=18, alpha=.6, color=c, label='%s (n=%d)' % (nm, len(d)))
ax[0].set_xlabel('Total-tumour / edema volume ratio'); ax[0].set_ylabel('Number of patients')
ax[0].set_xlim(0, 6); ax[0].legend(frameon=False, fontsize=9)
ax[0].set_title('(a) TT/ED ratio by MGMT status', fontweight='bold')
# (b) univariable AUC (natural direction)
cols = [r[0] for r in rows]; aus = [float(r[6]) for r in rows]
order = np.argsort(aus)
ax[1].barh([cols[i] for i in order], [aus[i] for i in order], color='#4C72B0')
ax[1].axvline(.5, color='k', ls='--', lw=.8)
ax[1].axvline(.687, color='#C44E52', ls=':', lw=1.4, label='Zahirovic 2025 (0.687)')
for i, kk in enumerate(order):
    ax[1].text(aus[kk] + .008, i, '%.3f' % aus[kk], va='center', fontsize=8)
ax[1].set_xlim(.3, .85); ax[1].set_xlabel('Univariate AUC (higher value → MGMT+)')
ax[1].legend(frameon=False, fontsize=8, loc='lower right')
ax[1].set_title('(b) Volumetric markers vs external benchmark', fontweight='bold')
# (c) by vendor
vs2 = list(sub.keys()); av = [sub[v]['auc'] for v in vs2]
bar = ax[2].bar(['All (n=%d)' % int(np.isfinite(ratio).sum())] + vs2,
                [au_r] + av, color=['#8172B3'] + ['#4C72B0', '#DD8452'][:len(vs2)])
for b, v in zip(bar, [au_r] + av):
    ax[2].text(b.get_x()+b.get_width()/2, v+.012, '%.3f' % v, ha='center', fontsize=9)
ax[2].axhline(.5, color='k', ls='--', lw=.8)
ax[2].axhline(.687, color='#C44E52', ls=':', lw=1.2)
ax[2].text(len(vs2)-0.35, .695, 'Zahirovic 2025 (0.687)', color='#C44E52', fontsize=8)
ax[2].set_ylim(.3, .8); ax[2].set_ylabel('AUC')
ax[2].set_title('(c) TT/ED ratio AUC by scanner vendor', fontweight='bold')
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig8_volume.png')); plt.close()
log('\nfig8 + table3 DONE (v3, polarity fixed, native-grid volumes)')
