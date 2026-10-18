# -*- coding: utf-8 -*-
r"""[LEGACY / NOT REPORTED IN THE PAPER] 简化 ComBat 调和 —— 早期探索性实现

*** 本文件不是稿件所报告 ComBat 结果的实现。***

稿件 §2.7(6) 报告的方法是 neuroCombat 0.2.12（经验贝叶斯收缩），每外折在训练折内拟合、
经 `neuroCombatFromTraining` 变换测试折。该实现位于：
    ../02_nested_cross_validation/paper1_nested_cv.py  第 [7] 节

本文件实现的是**无经验贝叶斯收缩**的简化 location/scale 版本（Johnson 2007），
属于早期探索性尝试，其输出（predictions_combat.csv / metrics_combat.json 的旧版）
未被稿件采用。保留仅用于溯源。
"""
r"""ComBat 厂商调和附加实验 (D:\Anaconda base 运行)

设计:
- 简化 ComBat (location/scale, 无EB收缩, Johnson 2007), 作为可 fit/transform 的变换器
- 在嵌套 CV 的每个外折内: 仅用训练折拟合调和参数, 应用到测试折 -> 无泄漏
- UIH n=1 无法估计批内方差, 剔除 (n=90: Siemens 57 / GE 33)
- 特征过滤: 缺失率>20% 与零方差 (不做 KW 过滤, 调和本就针对厂商差异)
- 模型: logistic / xgboost 嵌套5折 + LOVO(调和后)
- 另报告: 调和后厂商KW显著特征数(对比调和前) 验证调和效果
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

# 基础过滤
X = X[X.columns[X.isna().mean() <= 0.20]]
X = X[X.columns[X.var() > 1e-12]]
# 剔除 UIH
keep_m = vendor != 'UIH'
X, y, vendor = X.loc[keep_m], y[keep_m], vendor[keep_m]
print('X=%s y+=%d y-=%d vendors=%s' % (X.shape, y.sum(), (1-y).sum(),
      dict(pd.Series(vendor).value_counts())), flush=True)


class ComBat:
    """简化ComBat: fit(X_df, batch_arr) -> transform(X_df, batch_arr)"""
    def fit(self, Xdf, batch):
        self.feats = list(Xdf.columns)
        A = Xdf.values.astype(np.float64)
        self.batches = sorted(set(batch))
        n_all = len(A)
        # grand mean (按批加权) 与 pooled sd
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
            # 保留测试批自身的位置差异用于校准到训练批参数:
            # 用测试批相对训练批的均值偏移修正 gamma
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
    ranks = Xtr[order].rank()  # 预计算
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

# LOVO 调和后
lovo = {}
for v in ['Siemens', 'GE']:
    te_m = vendor == v
    p, _ = fit_predict('xgboost', X.iloc[~te_m], y[~te_m], X.iloc[te_m],
                       vendor[~te_m], vendor[te_m])
    lovo[v] = {'n_test': int(te_m.sum()), 'auc': float(roc_auc_score(y[te_m], p))}
    print('ComBat LOVO %s AUC=%.3f' % (v, lovo[v]['auc']), flush=True)
results['lovo_xgboost'] = lovo

# 调和效果验证: 调和前后厂商KW显著特征数
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
