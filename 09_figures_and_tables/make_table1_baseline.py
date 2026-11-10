# -*- coding: utf-8 -*-
r"""Paper 1 - Table 1 baseline table v2 (revised after peer review)

Changes from v1:
  1. Volumes now come from analysis/per_case_volumes.csv (native-grid voxel counting, n=89, no missing),
     no longer from PyRadiomics shape features (which held 2 spurious NaN cases).
  2. Table 1 gains the clinical variables the reviewers asked for: KPS, extent of resection,
     tumour location (lobe) and field strength; and a between-group p value for every variable
     (the reviewers called the absence of p values in v1 selective reporting).
  3. All variables now report their denominator.
  4. The study population is IDH wild-type GBM (WHO 2021 definition) throughout, so Table 1
     reports "IDH wild-type - 89 (100.0%)" without a between-group p value; TERT status and
     follow-up measures (overall survival / progression-free survival / recurrence) are outside
     the scope of Paper 1 and are not reported in Table 1 or used in any analysis.
  5. Output: table1_baseline.csv / .md (main table) and tableS5_baseline_full.csv (all variables, with p values).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id, resolve   # noqa: E402  patient code <-> local name
import os, json, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57'
CLIN = r'D:/1桌面文件/AI/GBM_clinical_data.xlsx'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
os.makedirs(ANA, exist_ok=True)
# Exclude the single United Imaging 1.5T case (see manuscript 3.8). Decided by vendor, not by name.
UIH_VENDOR = 'UIH'
SEED = 42

# ---------- cohort ----------
co = pd.read_csv(os.path.join(SRC, 'scripts', 'cohort_MGMT.csv'), encoding='utf-8-sig')
co = co[co.vendor.astype(str).str.strip() != UIH_VENDOR].copy()
co['patient'] = [to_id(p) for p in co['patient']]     # cohort table keys are Pxxx
# main analysis cohort: exactly the modelling cohort (cases with a usable feature matrix)
df = pd.read_csv(os.path.join(SRC, 'scripts', 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
df['patient'] = [to_id(p) for p in df['patient']]   # all data keys are Pxxx
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
X = X[X.columns[X.isna().mean() <= 0.20]]
X = X[X.columns[X.var() > 1e-12]]
X = X.loc[:, X.notna().any()]
common = [p for p in co.patient if p in X.index]
co = co.set_index('patient').loc[common].reset_index()

# ---------- clinical database ----------
cl = pd.read_excel(CLIN, sheet_name='GBM_clinical_data')
cl = cl.rename(columns={'Name': 'patient'})
cl['patient'] = cl['patient'].astype(str).str.strip()
cl['patient'] = [resolve(p) for p in cl['patient']]   # spelling variants in the clinical table -> Pxxx
cl = cl.drop_duplicates('patient', keep='last').set_index('patient')
M = co.merge(cl, on='patient', how='left', suffixes=('', '_cl'), validate='one_to_one')

# age and sex come from cohort_MGMT (so that they match the modelling definition)
M['age'] = pd.to_numeric(M['age'], errors='coerce')
M['male'] = M['gender'].astype(str).str.strip().str.upper().eq('M')
M['kps'] = pd.to_numeric(M['KPS'], errors='coerce')
M['eortc'] = M['Extent of resection'].astype(str).str.strip()
M['lobe'] = M['Location'].astype(str).str.strip()
M['mgmt_pos'] = M['MGMT'].astype(str).str.strip().eq('+')
M['side_l'] = M['site'].astype(str).str.strip().eq('Left')
M['side_r'] = M['site'].astype(str).str.strip().eq('Right')

# ---------- volumes: native-grid voxel counting ----------
pv = pd.read_csv(os.path.join(ANA, 'per_case_volumes.csv'))
pv['patient'] = pv['patient'].astype(str).str.strip()
pv = pv.drop_duplicates('patient').set_index('patient')
pv.index = [to_id(p) for p in pv.index]            # written tables are normalised to Pxxx
for k, c in [('tc', 'TC'), ('ne', 'NE'), ('ed', 'ED'), ('fl', 'FL'), ('br', 'BR')]:
    M[k] = M['patient'].map(pv[c])
M['tt'] = M['tc'] + M['ne']

# Pre-operative interval (surgery date - imaging date): taken from the interval already verified
# in the supplement (median 3 days, IQR 2-5). The value comes from the QC report; it is summarised here only, not recomputed.
PREOP_MED, PREOP_Q1, PREOP_Q3 = 3.0, 2.0, 5.0
PREOP_N = 89

# ---------- statistics helpers ----------
def mw(a, b):
    a = pd.Series(a).dropna().astype(float); b = pd.Series(b).dropna().astype(float)
    if len(a) < 3 or len(b) < 3: return np.nan
    return stats.mannwhitneyu(a, b, alternative='two-sided').pvalue

def cat_p(a, b):
    """a, b: 0/1 indicator Series (non-missing values only)"""
    a = pd.Series(a).dropna(); b = pd.Series(b).dropna()
    t = [[int(a.sum()), int(len(a) - a.sum())], [int(b.sum()), int(len(b) - b.sum())]]
    try:
        if min(min(t[0]), min(t[1])) < 5:
            return stats.fisher_exact(t).pvalue
        return stats.chi2_contingency(t, correction=False)[1]
    except Exception:
        return np.nan

def fmt_ms(s, d=1):
    s = pd.to_numeric(pd.Series(s), errors='coerce').dropna()
    if not len(s): return 'NA', 0
    return '%.*f ± %.*f' % (d, s.mean(), d, s.std()), len(s)

def fmt_med(s, d=1):
    s = pd.to_numeric(pd.Series(s), errors='coerce').dropna()
    if not len(s): return 'NA', 0
    q1, q3 = np.percentile(s, [25, 75])
    return '%.*f [%.*f–%.*f]' % (d, s.median(), d, q1, d, q3), len(s)

def fmt_cat(mask, denom_mask=None):
    m = pd.Series(mask).fillna(False).astype(bool)
    if denom_mask is None:
        den = m
    else:
        den = pd.Series(denom_mask).fillna(False).astype(bool)
    n_den = int(den.sum())
    if n_den == 0: return 'NA (0)'
    return '%d (%.1f%%)' % (int(m.sum()), 100 * m.sum() / n_den)

S = M['vendor'].astype(str).str.strip().eq('Siemens')
G = M['vendor'].astype(str).str.strip().eq('GE')

rows = []
def add(label, vals, denom=None, pval=np.nan, test=''):
    rows.append({'Characteristic': label, 'All (n=%d)' % len(M): vals[0],
                 'Siemens (n=%d)' % int(S.sum()): vals[1],
                 'GE (n=%d)' % int(G.sum()): vals[2],
                 'p (Siemens vs GE)': ('' if (isinstance(pval, float) and np.isnan(pval)) else '%.3f' % pval),
                 'Test': test, 'n reported': vals[3]})

# --- demographics ---
v, n = fmt_med(M['age']); rows.append({'Characteristic': 'Age, years (median [IQR])',
    'All (n=%d)' % len(M): v, 'Siemens (n=%d)' % int(S.sum()): fmt_med(M.loc[S, 'age'])[0],
    'GE (n=%d)' % int(G.sum()): fmt_med(M.loc[G, 'age'])[0],
    'p (Siemens vs GE)': '%.3f' % mw(M.loc[S, 'age'], M.loc[G, 'age']), 'Test': 'Mann–Whitney U', 'n reported': n})
rows.append({'Characteristic': 'Male sex — n (%)',
    'All (n=%d)' % len(M): fmt_cat(M['male'], M['gender'].notna()),
    'Siemens (n=%d)' % int(S.sum()): fmt_cat(M.loc[S, 'male'], M.loc[S, 'gender'].notna()),
    'GE (n=%d)' % int(G.sum()): fmt_cat(M.loc[G, 'male'], M.loc[G, 'gender'].notna()),
    'p (Siemens vs GE)': '%.3f' % cat_p(M.loc[S, 'male'].astype(float), M.loc[G, 'male'].astype(float)),
    'Test': 'χ²', 'n reported': int(M['gender'].notna().sum())})
v, n = fmt_med(M['kps'], 0); rows.append({'Characteristic': 'KPS (median [IQR])',
    'All (n=%d)' % len(M): v, 'Siemens (n=%d)' % int(S.sum()): fmt_med(M.loc[S, 'kps'], 0)[0],
    'GE (n=%d)' % int(G.sum()): fmt_med(M.loc[G, 'kps'], 0)[0],
    'p (Siemens vs GE)': '%.3f' % mw(M.loc[S, 'kps'], M.loc[G, 'kps']), 'Test': 'Mann–Whitney U', 'n reported': n})

# --- molecular markers (MGMT only; the whole cohort is IDH wild-type and TERT is outside the scope of this study) ---
den_mgmt = M['MGMT'].notna() & M['MGMT'].astype(str).str.strip().ne('nan')
p = cat_p(M.loc[S, 'mgmt_pos'].where(den_mgmt[S]).astype(float),
          M.loc[G, 'mgmt_pos'].where(den_mgmt[G]).astype(float))
rows.append({'Characteristic': 'MGMT promoter methylated — n (%)',
    'All (n=%d)' % len(M): fmt_cat(M['mgmt_pos'], den_mgmt),
    'Siemens (n=%d)' % int(S.sum()): fmt_cat(M.loc[S, 'mgmt_pos'], den_mgmt[S]),
    'GE (n=%d)' % int(G.sum()): fmt_cat(M.loc[G, 'mgmt_pos'], den_mgmt[G]),
    'p (Siemens vs GE)': ('%.3f' % p) if not np.isnan(p) else '',
    'Test': 'Fisher exact', 'n reported': int(den_mgmt.sum())})

# IDH: the whole cohort is IDH wild-type (WHO 2021 GBM definition) -> no between-group comparison
rows.append({'Characteristic': 'IDH wild-type (WHO 2021 GBM) — n (%)',
    'All (n=%d)' % len(M): '%d (100.0%%)' % len(M),
    'Siemens (n=%d)' % int(S.sum()): '%d (100.0%%)' % int(S.sum()),
    'GE (n=%d)' % int(G.sum()): '%d (100.0%%)' % int(G.sum()),
    'p (Siemens vs GE)': '—', 'Test': '—', 'n reported': len(M)})

# --- tumour side (all 89 cases: Left/Right/Bilateral/Midline) ---
SITE_ALL = ['Left', 'Right', 'Bilateral', 'Midline']
SDEN = M['site'].astype(str).str.strip().isin(SITE_ALL)
for lab, val in [('Tumour side — left — n (%)', 'Left'),
                 ('Tumour side — right — n (%)', 'Right'),
                 ('Tumour side — bilateral/midline — n (%)', '|')]:
    if val == '|':
        mask = M['site'].astype(str).str.strip().isin(['Bilateral', 'Midline'])
    else:
        mask = M['site'].astype(str).str.strip().eq(val)
    p = cat_p(mask[S].astype(float), mask[G].astype(float))
    rows.append({'Characteristic': lab,
        'All (n=%d)' % len(M): fmt_cat(mask, SDEN),
        'Siemens (n=%d)' % int(S.sum()): fmt_cat(mask[S], SDEN[S]),
        'GE (n=%d)' % int(G.sum()): fmt_cat(mask[G], SDEN[G]),
        'p (Siemens vs GE)': ('%.3f' % p) if not np.isnan(p) else '',
        'Test': 'χ²', 'n reported': int(SDEN.sum())})

# --- extent of resection ---
for k in ['GTR', 'PTR', 'BTR']:
    mask = M['eortc'].eq(k)
    den = M['eortc'].isin(['GTR', 'PTR', 'BTR'])
    p = cat_p(mask[S].astype(float), mask[G].astype(float))
    rows.append({'Characteristic': 'Extent of resection — %s — n (%%)' % k,
        'All (n=%d)' % len(M): fmt_cat(mask, den),
        'Siemens (n=%d)' % int(S.sum()): fmt_cat(mask[S], den[S]),
        'GE (n=%d)' % int(G.sum()): fmt_cat(mask[G], den[G]),
        'p (Siemens vs GE)': ('%.3f' % p) if not np.isnan(p) else '', 'Test': 'χ² / Fisher', 'n reported': int(den.sum())})

# --- tumour location (per lobe involvement, non-exclusive; revised 2026-11-04) ---
# The original version reported only "temporal lobe involvement" and "multilobar / corpus callosum
# (slash present)", which missed the frontal, parietal and occipital lobes altogether; and "multilobar"
# was really a proxy for "the report lists >= 2 sites separated by a slash", mixed into the corpus
# callosum row. It is now: four cortical lobes plus corpus callosum, one row each (counted whenever
# involved, so one case can appear under several lobes), and a separate "multilobar (>= 2 sites)" row
# to match the multi flag used in modelling.
den = M['lobe'].ne('nan') & M['lobe'].notna()
for _lab, _kw in [('frontal lobe', 'Frontal lobe'), ('temporal lobe', 'Temporal lobe'),
                  ('parietal lobe', 'Parietal lobe'), ('occipital lobe', 'Occipital lobe'),
                  ('corpus callosum', 'Corpus callosum')]:
    mask = M['lobe'].str.contains(_kw, case=False, na=False)
    p = cat_p(mask[S].astype(float), mask[G].astype(float))
    rows.append({'Characteristic': 'Tumour location — %s involvement — n (%%)' % _lab,
        'All (n=%d)' % len(M): fmt_cat(mask, den),
        'Siemens (n=%d)' % int(S.sum()): fmt_cat(mask[S], den[S]),
        'GE (n=%d)' % int(G.sum()): fmt_cat(mask[G], den[G]),
        'p (Siemens vs GE)': ('%.3f' % p) if not np.isnan(p) else '', 'Test': 'χ² / Fisher', 'n reported': int(den.sum())})
multi = M['lobe'].str.contains('/', na=False)
p = cat_p(multi[S].astype(float), multi[G].astype(float))
rows.append({'Characteristic': 'Tumour location — multilobar (≥2 sites reported) — n (%)',
    'All (n=%d)' % len(M): fmt_cat(multi, den),
    'Siemens (n=%d)' % int(S.sum()): fmt_cat(multi[S], den[S]),
    'GE (n=%d)' % int(G.sum()): fmt_cat(multi[G], den[G]),
    'p (Siemens vs GE)': ('%.3f' % p) if not np.isnan(p) else '', 'Test': 'χ² / Fisher', 'n reported': int(den.sum())})

# --- volumes ---
for lab, col in [('Tumour core volume, mL', 'tc'), ('Necrosis volume, mL', 'ne'),
                 ('Oedema (SNFH) volume, mL', 'ed'), ('FLAIR changes volume, mL', 'fl'),
                 ('Total tumour volume, mL', 'tt'), ('Brain volume, L', 'br')]:
    d = 2 if col == 'br' else 1
    ss = M[col] if col != 'br' else M['br'] / 1000.0
    v, n = fmt_med(ss, d)
    rows.append({'Characteristic': lab + ' (median [IQR])',
        'All (n=%d)' % len(M): v, 'Siemens (n=%d)' % int(S.sum()): fmt_med(ss[S], d)[0],
        'GE (n=%d)' % int(G.sum()): fmt_med(ss[G], d)[0],
        'p (Siemens vs GE)': '%.3f' % mw(ss[S], ss[G]), 'Test': 'Mann–Whitney U', 'n reported': n})

# --- acquisition ---
rows.append({'Characteristic': 'Field strength 3.0 T — n (%)',
    'All (n=%d)' % len(M): '%d (100.0%%)' % len(M),
    'Siemens (n=%d)' % int(S.sum()): '%d (100.0%%)' % int(S.sum()),
    'GE (n=%d)' % int(G.sum()): '%d (100.0%%)' % int(G.sum()),
    'p (Siemens vs GE)': '—', 'Test': '—', 'n reported': len(M)})
thin = M['patient'].isin(set())  # placeholder: slice thickness lives in Supplementary Table S1
rows.append({'Characteristic': 'T1c 3D thin-slice (1.0–1.2 mm) — n (%)',
    'All (n=%d)' % len(M): '59 (66.3%)',
    'Siemens (n=%d)' % int(S.sum()): '28 (49.1%)',
    'GE (n=%d)' % int(G.sum()): '31 (96.9%)',
    'p (Siemens vs GE)': '<0.001', 'Test': 'χ²', 'n reported': 89})
rows.append({'Characteristic': 'Preoperative MRI-to-surgery interval, days (median [IQR])',
    'All (n=%d)' % len(M): '%.0f [%.0f–%.0f]' % (PREOP_MED, PREOP_Q1, PREOP_Q3),
    'Siemens (n=%d)' % int(S.sum()): 'see Supplementary QC report', 'GE (n=%d)' % int(G.sum()): 'see Supplementary QC report',
    'p (Siemens vs GE)': '', 'Test': '—', 'n reported': PREOP_N})

t1 = pd.DataFrame(rows)

t1.to_csv(os.path.join(ANA, 'table1_baseline.csv'), index=False, encoding='utf-8-sig')
t1.drop(columns=['n reported']).to_csv(os.path.join(ANA, 'tableS5_baseline_full.csv'),
                                       index=False, encoding='utf-8-sig')
# main Table 1 (for the manuscript: Test column dropped, p values kept)
t1.drop(columns=['n reported', 'Test']).to_csv(os.path.join(ANA, 'table1_main.csv'),
                                               index=False, encoding='utf-8-sig')
t1['p (Siemens vs GE)'] = [('%s %s' % (r['p (Siemens vs GE)'], ('(%s)' % r['Test']) if r['Test'] not in ('—', '') else '')).strip()
                           for _, r in t1.iterrows()]
with open(os.path.join(ANA, 'table1_baseline.md'), 'w', encoding='utf-8') as f:
    f.write(t1.drop(columns=['n reported']).to_markdown(index=False))
    f.write('\n\n### Siemens vs GE tests (all variables)\n')
    for _, r in t1.iterrows():
        if r['p (Siemens vs GE)'] not in ('', '—'):
            f.write('- %s: p = %s (%s)\n' % (r['Characteristic'], r['p (Siemens vs GE)'], r['Test']))

print(t1.drop(columns=['Test']).to_string(index=False))
print('\nDONE table1 v2')

# ---------- Table 3 (condensed version for the manuscript) ----------
t3 = pd.read_csv(os.path.join(ANA, 'table3_volume.csv'))
ren = {'Marker': 'Marker', 'Higher in': 'Higher in', 'MGMT-unmethylated median [IQR]': 'MGMT− (n=44)',
       'MGMT-methylated median [IQR]': 'MGMT+ (n=45)', 'p (Mann–Whitney)': 'p',
       'AUC (predict MGMT+)': 'AUC', 'AUC 95% CI': '95% CI', 'BH-FDR q (7-marker family)': 'q'}
t3m = t3.rename(columns={k: v for k, v in ren.items() if k in t3.columns})
t3m = t3m.drop(columns=[c for c in ['n'] if c in t3m.columns])
if 'p' in t3m.columns: t3m['p'] = t3m['p'].map(lambda x: '%.3f' % x)
if 'q' in t3m.columns: t3m['q'] = t3m['q'].map(lambda x: '%.3f' % x)
if 'AUC' in t3m.columns: t3m['AUC'] = t3m['AUC'].map(lambda x: '%.3f' % x)
t3m.to_csv(os.path.join(ANA, 'table3_main.csv'), index=False, encoding='utf-8-sig')
print(t3m.to_string(index=False))

# key summary numbers for cross-checking the text
summ = {
    'n': int(len(M)), 'siemens': int(S.sum()), 'ge': int(G.sum()),
    'mgmt_pos': int(M['mgmt_pos'].sum()),
    'age_median': float(M['age'].median()), 'age_sd': float(M['age'].std()),
    'kps_median': float(M['kps'].median()),
    'tc_mean': float(M['tc'].mean()), 'ne_mean': float(M['ne'].mean()),
    'ed_mean': float(M['ed'].mean()), 'fl_mean': float(M['fl'].mean()),
    'tc_median': float(M['tc'].median()), 'ne_median': float(M['ne'].median()),
    'ed_median': float(M['ed'].median()), 'fl_median': float(M['fl'].median()),
    'eortc': {k: int(M['eortc'].eq(k).sum()) for k in ['GTR', 'PTR', 'BTR']},
    'idh_wildtype_n': int(len(M)),   # whole cohort IDH wild-type (WHO 2021)
    'kps_n': int(M['kps'].notna().sum()),
    'eortc_n': int(M['eortc'].isin(['GTR', 'PTR', 'BTR']).sum()),
    'lobe_n': int(M['lobe'].notna().sum()),
    'site_n': int(M['site'].astype(str).str.strip().isin(['Left', 'Right', 'Bilateral', 'Midline']).sum()),
}
json.dump(summ, open(os.path.join(ANA, 'table1_summary.json'), 'w'), indent=2)
print(json.dumps(summ, indent=2))
