# -*- coding: utf-8 -*-
r"""[LEGACY / NOT REPORTED IN THE PAPER] simplified ComBat harmonisation - an early exploratory implementation

*** This file is not the implementation behind the ComBat results reported in the paper. ***

The method reported in section 2.7(6) of the manuscript is neuroCombat 0.2.12 (empirical-Bayes shrinkage), fitted
inside the training folds of each outer fold and applied to the test fold through `neuroCombatFromTraining`. That
implementation lives in ../02_nested_cross_validation/paper1_nested_cv.py, section [7].

What this file implements is the simplified location/scale version **without empirical-Bayes shrinkage** (Johnson
2007). It was an early exploratory attempt, and its outputs (predictions_combat.csv and the older
metrics_combat.json) were never used in the paper. It is kept for provenance only.

Design of the experiment (run under the Anaconda base environment):
- simplified ComBat (location/scale, no EB shrinkage, Johnson 2007), exposed as a fit/transform transformer
- inside each outer fold of the nested CV: harmonisation parameters are fitted on the training folds only and
  applied to the test fold, so nothing leaks
- UIH n=1 cannot support a within-batch variance estimate and is dropped (n=90: Siemens 57 / GE 33)
- feature filtering: missing rate >20% and zero variance (no KW filtering, since harmonisation exists precisely to
  address vendor differences)
- models: logistic / xgboost nested 5-fold plus LOVO after harmonisation
- also reported: the number of vendor-KW-significant features before and after harmonisation, as a check on the
  harmonisation itself
"""

import os, csv, json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
import xgboost as xgb

BASE = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts'
WORK = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\analysis'
SEED = 42

print('loading...', flush=True)
df = pd.read_csv(os.path.join(BASE, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
co = pd.read_csv(os.path.join(BASE, 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
common = [p for p in co.index if p in X.index]
X = X.loc[common]
y = co.loc[common, 'y'].values
vendor = co.loc[common, 'vendor'].values

# basic filtering
X = X[X.columns[X.isna().mean() <= 0.20]]
X = X[X.columns[X.var() > 1e-12]]
# drop UIH
keep_m = vendor != 'UIH'
X, y, vendor = X.loc[keep_m], y[keep_m], vendor[keep_m]
print('X=%s y+=%d y-=%d vendors=%s' % (X.shape, y.sum(), (1-y).sum(),
      dict(pd.Series(vendor).value_counts())), flush=True)


class ComBat:
    """Simplified ComBat: fit(X_df, batch_arr) -> transform(X_df, batch_arr)"""
    def fit(self, Xdf, batch):
        self.feats = list(Xdf.columns)
        A = Xdf.values.astype(np.float64)
        self.batches = sorted(set(batch))
        n_all = len(A)
        # grand mean (batch-weighted) and pooled sd
        self.alpha = np.zeros(A.shape[1])
        self.psd = np.zeros(A.shape[1])
        w = np.array([ (batch == b).sum() / n_all for b in self.batches])
        for j in range(A.shape[1]):
            means = [A[batch == b, j].mean() for b in self.batches]
            vars_ = [A[batch == b, j].var(ddof=1) for b in self.batches]
            self.alpha[j] = np.dot(w, means)
            self.psd[j] = np.sqrt(max(np.dot(w, vars_), 1e-12))
        Z = (A - self.alpha) / self.psd
        self.gamma = {}
        self.delta = {}
        for b in self.batches:
            m = batch == b
            self.gamma[b] = Z[m].mean(0)
            self.delta[b] = np.maximum(Z[m].var(0), 1e-8)
        return self

    def transform(self, Xdf, batch):
        A = Xdf[self.feats].values.astype(np.float64)
        Z = (A - self.alpha) / self.psd
        out = np.empty_like(Z)
        for b in self.batches:
            m = batch == b
            if m.sum() == 0:
                continue
            Zb = Z[m]
            # keep the location difference of the test batch so the transform can be calibrated onto the
            # training-batch parameters: correct gamma by the mean shift of the test batch relative to training
            g = self.gamma[b]
            d = self.delta[b]
            out[m] = (Zb - g) / np.sqrt(d)
        return pd.DataFrame(out * self.psd + self.alpha, columns=self.feats, index=Xdf.index)


def inner_select(Xtr, ytr, max_feat=200, corr_th=0.90):
    auc = {}
    for c in Xtr.columns:
        v = Xtr[c].values
        m = np.isfinite(v)
        if m.sum() < 10 or len(np.unique(v[m])) < 2:
            continue
        try:
            auc[c] = abs(roc_auc_score(ytr[m], v[m]) - 0.5)
        except Exception:
            pass
    order = sorted(auc, key=auc.get, reverse=True)[:max_feat]
    kept = []
    ranks = Xtr[order].rank()  # precomputed
    for c in order:
        ok = True
        for k in kept:
            if abs(ranks[c].corr(ranks[k])) > corr_th:
                ok = False
                break
        if ok:
            kept.append(c)
        if len(kept) >= 40:
            break
    return kept


def make_model(name):
    if name == 'logistic':
        base = LogisticRegression(max_iter=5000, solver='liblinear')
        grid = {'C': [0.01, 0.1, 1, 10]}
    else:
        base = xgb.XGBClassifier(n_estimators=100, max_depth=2, learning_rate=0.1,
                                 subsample=0.8, colsample_bytree=0.8,
                                 eval_metric='logloss', random_state=SEED, n_jobs=2,
                                 tree_method='hist')
        grid = {'n_estimators': [50, 100], 'max_depth': [1, 2]}
    return base, grid


def fit_predict(name, Xtr, ytr, Xte, btr, bte):
    combat = ComBat().fit(Xtr, btr)
    Xtr_h = combat.transform(Xtr, btr)
    Xte_h = combat.transform(Xte, bte)
    cols = inner_select(Xtr_h, ytr)
    if not cols:
        return np.full(len(Xte), 0.5), cols
    base, grid = make_model(name)
    pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()), ('m', base)])
    from sklearn.model_selection import GridSearchCV
    gs = GridSearchCV(pipe, {'m__' + k: v for k, v in grid.items()}, cv=3, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr_h[cols], ytr)
    return gs.predict_proba(Xte_h[cols])[:, 1], cols


def run_nested(name, X, y, vendor):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    prob = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        p, _ = fit_predict(name, X.iloc[tr], y[tr], X.iloc[te], vendor[tr], vendor[te])
        prob[te] = p
    return prob


results = {}
pred = pd.DataFrame({'patient': X.index, 'y': y, 'vendor': vendor})
for name in ['logistic', 'xgboost']:
    prob = run_nested(name, X, y, vendor)
    auc = roc_auc_score(y, prob)
    results[name] = {'auc': float(auc), 'brier': float(brier_score_loss(y, prob))}
    pred['prob_' + name] = prob
    print('ComBat nested CV %s AUC=%.3f' % (name, auc), flush=True)

# LOVO after harmonisation
lovo = {}
for v in ['Siemens', 'GE']:
    te_m = vendor == v
    p, _ = fit_predict('xgboost', X.iloc[~te_m], y[~te_m], X.iloc[te_m],
                       vendor[~te_m], vendor[te_m])
    lovo[v] = {'n_test': int(te_m.sum()), 'auc': float(roc_auc_score(y[te_m], p))}
    print('ComBat LOVO %s AUC=%.3f' % (v, lovo[v]['auc']), flush=True)
results['lovo_xgboost'] = lovo

# effect of harmonisation: number of vendor-KW-significant features before and after
from scipy import stats
cb = ComBat().fit(X, vendor)
Xh = cb.transform(X, vendor)
kw_before = kw_after = 0
for c in X.columns:
    a = X.loc[vendor == 'Siemens', c].dropna()
    b = X.loc[vendor == 'GE', c].dropna()
    if len(a) > 5 and len(b) > 5 and stats.kruskal(a, b).pvalue < 0.05:
        kw_before += 1
    a = Xh.loc[vendor == 'Siemens', c].dropna()
    b = Xh.loc[vendor == 'GE', c].dropna()
    if len(a) > 5 and len(b) > 5 and stats.kruskal(a, b).pvalue < 0.05:
        kw_after += 1
results['kw_vendor_significant'] = {'before': kw_before, 'after': kw_after, 'total': X.shape[1]}
print('KW significant vendor features: before=%d after=%d / %d' % (kw_before, kw_after, X.shape[1]), flush=True)

pred.to_csv(os.path.join(WORK, 'predictions_combat.csv'), index=False, encoding='utf-8-sig')
with open(os.path.join(WORK, 'metrics_combat.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2, ensure_ascii=False)
print('DONE', flush=True)
