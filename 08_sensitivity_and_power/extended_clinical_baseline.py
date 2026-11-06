# -*- coding: utf-8 -*-
r"""Paper 1 -- R2-Major5 (a stronger clinical baseline) + R1-M6 (calibration slope and intercept)

A. Full clinical baseline: age + sex + temporal-lobe involvement + tumour volume (TumorCore) [+ oedema volume]
   (no feature selection, straight logistic regression; equivalent to the repeated 5x5 stratified nested CV,
   implemented as repeated 5-fold x5)
B. Calibration: fit logit(p) ~ logit(y) on the pooled OOF predictions of the radiomics logistic model and report
   the calibration slope and intercept (intercept-in-the-large).

Output: analysis/extra_clinical_calibration.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id, resolve   # noqa: E402  patient code <-> local name
import os, json, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, brier_score_loss
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
CLIN = r'D:/1桌面文件/AI/GBM_clinical_data.xlsx'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
SEED = 42
REPEATS = 5

co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig')
# Exclude the single United Imaging 1.5T case (see manuscript 3.8). Decided by vendor, not by name:
# in the released data patients are Pxxx, so excluding by name fails silently and changes n.
co = co[co.vendor.astype(str).str.strip() != 'UIH'].set_index('patient')
co.index = [to_id(p) for p in co.index]            # all data keys are Pxxx
co['y'] = co['MGMT'].astype(str).str.strip().str.startswith('+').astype(int)

pv = pd.read_csv(os.path.join(ANA, 'per_case_volumes.csv'))
pv['patient'] = pv['patient'].astype(str).str.strip()
pv = pv.drop_duplicates('patient').set_index('patient')
pv.index = [to_id(p) for p in pv.index]            # written tables are normalised to Pxxx

cl = pd.read_excel(CLIN, sheet_name='GBM_clinical_data').rename(columns={'Name': 'patient'})
# the Name column may carry a parenthetical note and a few cases are spelt differently from the imaging directory -> normalise to codes with resolve() first
cl['patient'] = (cl['patient'].astype(str).str.replace(r'\(.*?\)', '', regex=True).str.strip())
cl = cl.drop_duplicates('patient').set_index('patient')
cl.index = [resolve(p) for p in cl.index]          # same namespace as co / pv (Pxxx)

idx = [p for p in co.index if p in pv.index and p in cl.index]
D = pd.DataFrame(index=idx)
D['y'] = co.loc[idx, 'y']
D['age'] = pd.to_numeric(co.loc[idx, 'age'], errors='coerce')
D['male'] = co.loc[idx, 'gender'].astype(str).str.strip().str.upper().eq('M').astype(float)
D['temporal'] = cl.loc[idx, 'Location'].astype(str).str.contains('Temporal', case=False, na=False).astype(float)
D['tcv'] = pv.loc[idx, 'TC'].astype(float).values
D['edv'] = pv.loc[idx, 'ED'].astype(float).values
D['multilobar'] = cl.loc[idx, 'Location'].astype(str).str.contains('/', na=False).astype(float)
print('n =', len(D), '| MGMT+ =', int(D.y.sum()))

def rep_cv_auc(cols, D, repeats=REPEATS, seed=SEED):
    aucs = []
    for r in range(repeats):
        skf = StratifiedKFold(5, shuffle=True, random_state=seed + r)
        p = np.zeros(len(D))
        for tr, te in skf.split(D[cols], D['y']):
            pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()),
                             ('m', LogisticRegression(max_iter=5000))])
            pipe.fit(D[cols].iloc[tr], D['y'].iloc[tr])
            p[te] = pipe.predict_proba(D[cols].iloc[te])[:, 1]
        aucs.append(roc_auc_score(D['y'], p))
    return float(np.mean(aucs)), float(np.std(aucs)), aucs

res = {}
for name, cols in [
    ('age+sex', ['age', 'male']),
    ('age+sex+location(temporal/multilobar)', ['age', 'male', 'temporal', 'multilobar']),
    ('age+sex+tumourcore_volume', ['age', 'male', 'tcv']),
    ('age+sex+location+tumourcore_volume', ['age', 'male', 'temporal', 'multilobar', 'tcv']),
    ('age+sex+location+tumourcore+oedema_volume', ['age', 'male', 'temporal', 'multilobar', 'tcv', 'edv']),
]:
    m, s, a = rep_cv_auc(cols, D)
    res[name] = {'auc_mean': round(m, 4), 'auc_sd': round(s, 4), 'n_vars': len(cols)}
    print('%-45s AUC = %.4f ± %.4f (k=%d)' % (name, m, s, len(cols)))

# ---- B. calibration (radiomics logistic OOF) ----
pr = pd.read_csv(os.path.join(ANA, 'predictions.csv'))
print('\npredictions.csv columns:', pr.columns.tolist(), '| shape', pr.shape)
cols_low = [c.lower() for c in pr.columns]
cand = [c for c in pr.columns if 'logist' in c.lower() or c.lower() in ('p_logistic', 'logistic')]
if cand:
    pc = cand[0]
    if 'patient' in pr.columns:
        pr = pr.set_index('patient')
    pr.index = [to_id(p) for p in pr.index]        # written tables are normalised to Pxxx
    m = pr.index.isin(D.index)
    yy = D.loc[[p for p in pr.index[m]], 'y'].values
    pp = pr.loc[pr.index[m], pc].values.astype(float)
    ok = np.isfinite(pp) & np.isfinite(yy)
    yy, pp = yy[ok], pp[ok]
    eps = 1e-6
    lg = np.log(np.clip(pp, eps, 1 - eps) / (1 - np.clip(pp, eps, 1 - eps)))
    Xc = np.column_stack([np.ones(len(lg)), lg])
    lr = LogisticRegression(penalty=None, max_iter=5000, solver='lbfgs')
    lr.fit(Xc, yy)
    slope = float(lr.coef_[0][1]); intercept = float(lr.intercept_[0])
    res['calibration_radiomics_logistic_oof'] = {
        'n': int(len(yy)), 'calibration_slope': slope, 'calibration_intercept': intercept,
        'brier': float(brier_score_loss(yy, pp)), 'auc': float(roc_auc_score(yy, pp)),
        'column_used': pc, 'note': 'logit(p) ~ logit(y): slope 1 is perfect (0.1 is poor), intercept 0 is perfect'}
    print('calibration slope=%.3f intercept=%.3f n=%d' % (slope, intercept, len(yy)))
else:
    print('WARNING: no logistic probability column found')

json.dump(res, open(os.path.join(ANA, 'extra_clinical_calibration.json'), 'w'), indent=2)
print('\nsaved -> analysis/extra_clinical_calibration.json')
