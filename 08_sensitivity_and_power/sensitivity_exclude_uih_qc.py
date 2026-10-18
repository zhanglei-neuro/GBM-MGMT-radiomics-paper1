# -*- coding: utf-8 -*-
r"""
Paper 1 -- 敏感性分析（rev.7: 影像日期异常已修复, 队列口径更新）
1) 排除 United Imaging 病例 (P021, n=1)  --- 用户问题 1
2) 仅保留 3D 薄层 (<=1.2mm) 子集 --- 用户问题 2 的技术对照

完全复用 paper1_analysis.py 的建模管线 (同一 SEED / CV 结构 / 特征筛选)。
输出: analysis/sensitivity_cohorts.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import os, json, warnings
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
import xgboost as xgb
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
SEED, REPEATS, NSPLIT = 42, 5, 5

df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
Xall = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]        # 数据键统一为 Pxxx
Xall.index = [to_id(p) for p in Xall.index]
common = [p for p in co.index if p in Xall.index]
X = Xall.loc[common].copy()
y = co.loc[common, 'y'].values
vendor = co.loc[common, 'vendor'].values
age = pd.to_numeric(co.loc[common, 'age'], errors='coerce').values
male = co.loc[common, 'gender'].astype(str).str.strip().eq('M').astype(int).values


def inner_select(Xtr, ytr, max_feat=200, corr_th=0.90, max_keep=40):
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
    kept, rk = [], Xtr[order].rank()
    for c in order:
        if all(abs(rk[c].corr(rk[k])) <= corr_th for k in kept):
            kept.append(c)
        if len(kept) >= max_keep:
            break
    return kept


def make_model(name, seed=SEED):
    if name == 'logistic':
        return (LogisticRegression(max_iter=5000, solver='liblinear'), {'C': [0.01, 0.1, 1, 10]})
    if name == 'xgboost':
        return (xgb.XGBClassifier(n_estimators=100, max_depth=2, learning_rate=0.1,
                                  subsample=0.8, colsample_bytree=0.8, eval_metric='logloss',
                                  random_state=seed, n_jobs=2, tree_method='hist'),
                {'n_estimators': [50, 100], 'max_depth': [1, 2]})
    return (RandomForestClassifier(n_estimators=200, random_state=seed, n_jobs=2),
            {'max_depth': [2, 4, None], 'min_samples_leaf': [2, 5]})


def fit_predict(name, Xtr, ytr, Xte, seed=SEED):
    cols = inner_select(Xtr, ytr)
    if not cols:
        return np.full(len(Xte), 0.5), cols
    base, grid = make_model(name, seed)
    pipe = Pipeline([('imp', SimpleImputer(strategy='median')),
                     ('sc', StandardScaler()), ('m', base)])
    gs = GridSearchCV(pipe, {'m__' + k: v for k, v in grid.items()}, cv=3,
                      scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr[cols], ytr)
    return gs.predict_proba(Xte[cols])[:, 1], cols


def rep_nested(name, X_, y_, repeats=REPEATS, seed=SEED):
    aucs = []
    for r in range(repeats):
        skf = StratifiedKFold(n_splits=NSPLIT, shuffle=True, random_state=seed + r)
        p = np.zeros(len(y_))
        for tr, te in skf.split(X_, y_):
            pp, _ = fit_predict(name, X_.iloc[tr], y_[tr], X_.iloc[te], seed + r)
            p[te] = pp
        aucs.append(roc_auc_score(y_, p))
    return float(np.mean(aucs)), float(np.std(aucs))


def boot_ci(y_, p_, n=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    a = []
    for _ in range(n):
        idx = rng.integers(0, len(y_), len(y_))
        if len(np.unique(y_[idx])) < 2:
            continue
        a.append(roc_auc_score(y_[idx], p_[idx]))
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))


def run_cohort(mask, label):
    Xc = X.loc[mask].copy()
    yc = y[mask]
    miss = Xc.isna().mean()
    Xc = Xc[miss[miss <= 0.20].index]
    vv = Xc.var()
    Xc = Xc[vv[vv > 1e-12].index]
    n_feat = Xc.shape[1]
    res = {'label': label, 'n': int(len(yc)), 'n_pos': int(yc.sum()),
           'n_features': int(n_feat), 'models': {}}
    for name in ['logistic', 'xgboost', 'randomforest']:
        m, s = rep_nested(name, Xc, yc)
        res['models'][name] = {'auc': round(m, 4), 'sd': round(s, 4)}
    # 临床基线
    Xcl = pd.DataFrame({'age': age[mask], 'male': male[mask]}, index=Xc.index)
    p = np.zeros(len(yc))
    for r in range(REPEATS):
        skf = StratifiedKFold(n_splits=NSPLIT, shuffle=True, random_state=SEED + r)
        pa = np.zeros(len(yc))
        for tr, te in skf.split(Xcl, yc):
            pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()),
                             ('m', LogisticRegression(max_iter=2000))])
            pipe.fit(Xcl.iloc[tr], yc[tr])
            pa[te] = pipe.predict_proba(Xcl.iloc[te])[:, 1]
        p = pa
    lo, hi = boot_ci(yc, p)
    res['models']['clinical'] = {'auc': round(float(roc_auc_score(yc, p)), 4),
                                 'ci': [round(lo, 4), round(hi, 4)],
                                 'brier': round(float(brier_score_loss(yc, p)), 4)}
    print('[%s] n=%d  feas=%d' % (label, res['n'], n_feat), flush=True)
    for k, v in res['models'].items():
        print('    %-13s AUC=%s' % (k, v['auc']), flush=True)
    return res


# ---------- 定义各组 ----------
idx = pd.Index(X.index)
uiH = np.array([str(co.loc[p, 'vendor']).strip() == 'UIH'
                for p in idx])   # 按 vendor 判断, 不按姓名

# 3D 薄层子集
s1b = pd.read_csv(os.path.join(OUT, 'supplementary/tables/S1b_paper1_per_case.csv'))
s1b['Patient'] = [to_id(p) for p in s1b['Patient']]   # 与 idx 同一命名空间
t1c = s1b[s1b['Sequence role'] == 'T1C'].copy()
t1c['th'] = pd.to_numeric(t1c['Slice thickness (mm)'], errors='coerce')
thin = set(t1c.loc[t1c.th <= 1.2, 'Patient'])
thick = set(t1c.loc[t1c.th >= 3, 'Patient'])
is_thin = idx.isin(list(thin))
is_thick = idx.isin(list(thick))

# 人口学平衡
print('=== 2D 厚层 vs 3D 薄层 人口学 ===')
for nm, m in [('3D thin', is_thin), ('2D thick', is_thick)]:
    print('  %-9s n=%2d  age=%.1f  male=%.0f%%  MGMT+=%d' % (
        nm, m.sum(), np.nanmean(age[m]), 100 * np.nanmean(male[m]), int(y[m].sum())))
print('  age p=%.3f  sex p=%.3f' % (
    __import__('scipy.stats', fromlist=['x']).mannwhitneyu(age[is_thin], age[is_thick]).pvalue,
    __import__('scipy.stats', fromlist=['x']).chi2_contingency(
        [[int(male[is_thin].sum()), int((1 - male[is_thin]).sum())],
         [int(male[is_thick].sum()), int((1 - male[is_thick]).sum())]])[1]))

cohorts = [
    (np.ones(len(y), bool), 'A_full_incl_UIH_%d' % len(y)),
    (~uiH, 'B_primary_no_UIH_%d' % int((~uiH).sum())),
    (is_thin & ~uiH, 'C_thin3D_only_%d' % int((is_thin & ~uiH).sum())),
    (is_thick & ~uiH, 'D_thick2D_only_%d' % int((is_thick & ~uiH).sum())),
]

out = []
for m, lab in cohorts:
    r = run_cohort(m, lab)
    r['tag'] = lab[0]        # 与 merge_sensitivity_results.py 同一口径 (A/B/C/D)
    out.append(r)

json.dump(out, open(os.path.join(ANA, 'sensitivity_cohorts.json'), 'w'),
          indent=2, ensure_ascii=False)
print('\nsaved -> analysis/sensitivity_cohorts.json')
