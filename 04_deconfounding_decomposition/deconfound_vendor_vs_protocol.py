# -*- coding: utf-8 -*-
r"""Paper 1 -- R1-M4(i): what remains of the "vendor effect" once acquisition parameters are adjusted for

Main analysis (the registered definition, identical to §3.x of the manuscript, Figure S4 and R1-M4 of the response letter)

For every radiomic feature f, with the acquisition parameters of its source sequence as covariates, fit

    z(f) ~ vendor + 2D/3D + slice thickness + TR + TE + TI

take the p value of the t test on the vendor coefficient, compare it with the Kruskal-Wallis p value of
"vendor alone", and count how the number of BH-FDR q<0.05 features changes. That quantifies what is left
of the vendor effect after the acquisition parameters are removed.
Expected results (the values quoted in the manuscript): 3382 (81.5%) -> 259 (6.2%); paired 224 / 3158 (93.4%);
median partial R2 vendor|acquisition 0.016 vs acquisition|vendor 0.098 (about six-fold).

Sensitivity analysis (not the manuscript's headline number; archived in sensitivity_inplane_covariate only)

Adds **in-plane pixel spacing** to the covariate set above
(= (Pixel spacing row + Pixel spacing col)/2, taken from the per-case DICOM measurement of the sequence
the feature actually belongs to, Supplementary Table S1b). Why include it: in-plane spacing in this cohort
runs from 0.34 to 1.30 mm and correlates with vendor (r=0.375), so leaving it out could bias the estimated
net vendor effect. Under this definition adjusted q<0.05 = 362 (8.7%) and the direction of the conclusion is
unchanged (about nine-tenths of the vendor effect is still explained by the acquisition parameters), so the
manuscript keeps reporting the five-covariate registered definition.

Variance decomposition: the median partial R2 across all features for (vendor | acquisition) and
(acquisition | vendor) under the main definition, and for (in-plane | the rest) under the sensitivity one.

Output: analysis/m4_vendor_vs_protocol.json
"""
import os, json, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
TAB = os.path.join(OUT, 'supplementary', 'tables')
SEED = 42
rng = np.random.default_rng(SEED)

COV_BASE = ['type2D', 'thk', 'TR', 'TE', 'TI']          # registered definition (5 covariates)
COV_SENS = COV_BASE + ['inplane']                       # sensitivity definition (adds in-plane pixel spacing)


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


# ---------- feature matrix ----------
df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
X = X[X.columns[X.isna().mean() <= 0.20]]
X = X[X.columns[X.var() > 1e-12]]
X = X.loc[:, X.notna().any()]

co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig').set_index('patient')
# Exclude the single United Imaging case (n=1): no platform-level inference is possible, and a
# one-case ComBat batch cannot estimate batch parameters. Excluded by vendor rather than by name so that
# pseudonymisation (<PATIENT>) cannot silently break it and push n to 90 instead of 89.
common = [p for p in co.index if p in X.index and str(co.loc[p, 'vendor']).strip() != 'UIH']
X = X.loc[common]
y_vendor = (co.loc[common, 'vendor'] == 'Siemens').astype(int)   # 1 = Siemens
print('cohort after excluding UIH: n=%d' % len(common))

# ---------- per-case acquisition parameters (the sequence actually used for extraction) ----------
acq = pd.read_csv(os.path.join(TAB, 'S1b_paper1_per_case.csv'), encoding='utf-8-sig')
acq.columns = [c.strip() for c in acq.columns]
acq = acq.rename(columns={'Sequence role': 'seq', 'Acquisition type': 'atype'})
acq['Patient'] = acq['Patient'].astype(str).str.strip()
VCOLS = ['TR (ms)', 'TE (ms)', 'TI (ms)', 'Slice thickness (mm)',
         'Pixel spacing row (mm)', 'Pixel spacing col (mm)']
for c in VCOLS:
    acq[c] = pd.to_numeric(acq[c], errors='coerce')
acq['type2D'] = acq['atype'].astype(str).str.contains('2D').astype(float)
# in-plane pixel spacing = mean of the row and column spacing
acq['inplane'] = acq[['Pixel spacing row (mm)', 'Pixel spacing col (mm)']].mean(axis=1)

# one row per (patient, seq)
A = acq.set_index(['Patient', 'seq'])
print('acquisition rows:', len(acq), 'patients:', acq.Patient.nunique())
print(acq.groupby('seq')[VCOLS + ['type2D']].median().to_string())

# ---------- design matrix (parameters swapped per feature) ----------
SEQ_OF = {}
for c in X.columns:
    parts = c.split('|')
    s = parts[1]
    SEQ_OF[c] = 'T1C' if s == 'T1CE' else s      # features say T1CE, the acquisition table says T1C


def design_for(seq, idx):
    """Return the covariate matrix (n x k) for that sequence"""
    sub = A.xs(seq, level='seq') if seq in set(acq['seq']) else None
    if sub is None or len(sub) == 0:
        return None
    sub = sub.reindex(idx)
    D = pd.DataFrame(index=idx)
    D['type2D'] = sub['type2D'].values
    D['thk'] = sub['Slice thickness (mm)'].values
    D['TR'] = sub['TR (ms)'].values
    D['TE'] = sub['TE (ms)'].values
    D['TI'] = sub['TI (ms)'].values
    D['inplane'] = sub['inplane'].values
    return D


# check that the sequence names match
avail_roles = set(acq['seq'].astype(str))
need_roles = set(SEQ_OF.values())
print('feature seq roles:', sorted(need_roles), '| acquisition roles:', sorted(avail_roles))
missing_roles = need_roles - avail_roles
if missing_roles:
    print('WARNING missing roles:', missing_roles)


# standardise the continuous covariates
def zs(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd > 0 else 1.0)


idx = X.index
cov_cache = {}
for seq in need_roles:
    D = design_for(seq, idx)
    if D is None:
        continue
    for c in ['thk', 'TR', 'TE', 'TI', 'inplane']:
        D[c] = zs(D[c])
    cov_cache[seq] = D

Xv = X.values.astype(float)


def bh_run(covs):
    """For a given covariate set, return (adj_p, pr2_vendor, pr2_acq, pr2_inplane).
    The pr2_* values are only meaningful when the corresponding columns are in the model; in-plane returns all NaN when it is not among covs."""
    adj_p = np.full(X.shape[1], np.nan)
    pr2_vendor = np.full(X.shape[1], np.nan)   # vendor | acquisition
    pr2_acq = np.full(X.shape[1], np.nan)      # acquisition | vendor
    pr2_inplane = np.full(X.shape[1], np.nan)  # in-plane | vendor + other acquisition parameters
    has_ip = 'inplane' in covs

    def arsq(res, res_r):
        if res.rsquared is None or res_r.rsquared is None:
            return np.nan
        return (res.rsquared - res_r.rsquared) / (1 - res_r.rsquared) \
            if 1 - res_r.rsquared > 1e-9 else np.nan

    for j, col in enumerate(X.columns):
        seq = SEQ_OF[col]
        D = cov_cache.get(seq)
        if D is None:
            continue
        v = Xv[:, j]
        m = np.isfinite(v) & np.isfinite(D[covs].values).all(axis=1)
        if m.sum() < 30:
            continue
        try:
            yy = zs(pd.Series(v[m]))
            Dm = D.loc[m]
            base = Dm[covs].copy()
            full = base.copy(); full['vendor'] = y_vendor.loc[m].values
            rf = sm.OLS(yy.values, sm.add_constant(full.astype(float))).fit()
            rr = sm.OLS(yy.values, sm.add_constant(base.astype(float))).fit()
            adj_p[j] = rf.pvalues['vendor']
            pr2_vendor[j] = arsq(rf, rr)
            # acquisition | vendor
            rr2 = sm.OLS(yy.values,
                         sm.add_constant(pd.DataFrame(
                             {'vendor': y_vendor.loc[m].values}).astype(float))).fit()
            pr2_acq[j] = arsq(rf, rr2)
            if has_ip:
                ri = sm.OLS(yy.values,
                            sm.add_constant(full.drop(columns=['inplane']).astype(float))).fit()
                pr2_inplane[j] = arsq(rf, ri)
        except Exception:
            pass
    return adj_p, pr2_vendor, pr2_acq, pr2_inplane


# ---------- vendor-only KW p values (baseline) ----------
kw_p = np.full(X.shape[1], np.nan)
for j in range(X.shape[1]):
    v = Xv[:, j]; m = np.isfinite(v)
    if m.sum() < 10:
        continue
    a = v[m & (y_vendor.values == 1)]; b = v[m & (y_vendor.values == 0)]
    if len(a) < 5 or len(b) < 5:
        continue
    kw_p[j] = stats.kruskal(a, b).pvalue

kw_q = bh(np.nan_to_num(kw_p, nan=1.0))

results = {'n': int(len(idx)), 'n_features': int(X.shape[1]),
           'note': ('primary model (registered, matches manuscript & Figure S4): '
                    'z(feature) ~ vendor + type2D + thickness + TR + TE + TI')}

# ---------- main analysis: five-covariate registered definition ----------
adj_p, pr2_vendor, pr2_acq, _ = bh_run(COV_BASE)
adj_q = bh(np.nan_to_num(adj_p, nan=1.0))
valid = np.isfinite(kw_p) & np.isfinite(adj_p)

results['unadjusted_vendor_only'] = {
    'raw_p05': int(np.nansum(kw_p < 0.05)),
    'q05': int(np.nansum(kw_q < 0.05)),
    'pct_q05': float(100 * np.nansum(kw_q < 0.05) / len(kw_q)),
}
results['adjusted_for_acquisition'] = {
    'raw_p05': int(np.nansum(adj_p < 0.05)),
    'q05': int(np.nansum(adj_q < 0.05)),
    'pct_q05': float(100 * np.nansum(adj_q < 0.05) / len(adj_q)),
}
results['paired_on_valid'] = {
    'n_valid_features': int(valid.sum()),
    'q05_both': int(np.sum((kw_q[valid] < 0.05) & (adj_q[valid] < 0.05))),
    'q05_unadj_only': int(np.sum((kw_q[valid] < 0.05) & (adj_q[valid] >= 0.05))),
    'q05_adj_only': int(np.sum((kw_q[valid] >= 0.05) & (adj_q[valid] < 0.05))),
}
results['partial_r2'] = {
    'vendor_given_acquisition_median': float(np.nanmedian(pr2_vendor)),
    'acquisition_given_vendor_median': float(np.nanmedian(pr2_acq)),
    'vendor_given_acquisition_mean': float(np.nanmean(pr2_vendor)),
    'acquisition_given_vendor_mean': float(np.nanmean(pr2_acq)),
}

# ---------- sensitivity: add in-plane pixel spacing (six covariates) ----------
adj_p6, pr2_vendor6, pr2_acq6, pr2_inplane = bh_run(COV_SENS)
adj_q6 = bh(np.nan_to_num(adj_p6, nan=1.0))
valid6 = np.isfinite(kw_p) & np.isfinite(adj_p6)
results['sensitivity_inplane_covariate'] = {
    'note': ('sensitivity only, NOT the manuscript headline number: '
             'z(feature) ~ vendor + type2D + thickness + TR + TE + TI '
             '+ in-plane pixel spacing'),
    'raw_p05': int(np.nansum(adj_p6 < 0.05)),
    'q05': int(np.nansum(adj_q6 < 0.05)),
    'pct_q05': float(100 * np.nansum(adj_q6 < 0.05) / len(adj_q6)),
    'paired_on_valid': {
        'n_valid_features': int(valid6.sum()),
        'q05_both': int(np.sum((kw_q[valid6] < 0.05) & (adj_q6[valid6] < 0.05))),
        'q05_unadj_only': int(np.sum((kw_q[valid6] < 0.05) & (adj_q6[valid6] >= 0.05))),
        'q05_adj_only': int(np.sum((kw_q[valid6] >= 0.05) & (adj_q6[valid6] < 0.05))),
    },
    'partial_r2': {
        'vendor_given_acquisition_median': float(np.nanmedian(pr2_vendor6)),
        'acquisition_given_vendor_median': float(np.nanmedian(pr2_acq6)),
        'inplane_given_rest_median': float(np.nanmedian(pr2_inplane)),
        'inplane_given_rest_mean': float(np.nanmean(pr2_inplane)),
    },
}

# Collinearity diagnostics (using the T1C design matrix as the example; the pairwise Pearson correlations
# do not depend on whether in-plane is included). Note that SEQ_OF already maps the feature-side T1CE to the
# acquisition-table side T1C, so the keys of cov_cache are 'T1C'.
D = cov_cache['T1C'] if 'T1C' in cov_cache else cov_cache.get('T1CE')
if D is not None:
    full = D.copy(); full['vendor'] = y_vendor.values
    c = full.corr().abs()
    results['collinearity_max_abs_r_with_vendor'] = \
        {k: float(v) for k, v in c['vendor'].drop('vendor').items()}

json.dump(results, open(os.path.join(ANA, 'm4_vendor_vs_protocol.json'), 'w'), indent=2)
print(json.dumps(results, indent=2))
print('saved -> analysis/m4_vendor_vs_protocol.json')
