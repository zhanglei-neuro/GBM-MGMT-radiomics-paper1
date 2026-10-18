# -*- coding: utf-8 -*-
r"""Paper 1 -- R2-Major3: GE 亚组内 AUC 是否可由随机波动解释
对 GE 亚组 (n=32, MGMT+ 13 / - 19) 做标签置换 null 分布:
  - 置换标签 100 次, 每次用与主分析完全相同的管线 (单变量|AUC-.5| top200 -> 去相关
    -> <=40 特征 -> L2 logistic 内层 3 折调参) 跑 5 折 CV, 计算合并 OOF AUC
  - 报告 null 分布分位数, 并给出实际 AUC (logistic 0.328 / XGBoost 0.220) 的经验 p
输出: analysis/ge_perm_null.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import os, sys, json, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
SEED = 42
N_PERM = 100
def log(*a): print(*a, flush=True)

df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
Xall = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
Xall.index = [to_id(p) for p in Xall.index]        # 数据键统一为 Pxxx
co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]
# 排除 United Imaging 1.5T 单例 (见稿件 3.8)。按 vendor 判断, 不按姓名:
# 发布版数据的患者标识为 Pxxx, 按姓名排除会静默失效并改变 n。
common = [p for p in co.index if p in Xall.index
          and str(co.loc[p, 'vendor']).strip() != 'UIH']
X = Xall.loc[common].copy()
X = X[X.columns[X.isna().mean() <= 0.20]]; X = X[X.columns[X.var() > 1e-12]]
y = co.loc[X.index, 'y'].values
vendor = co.loc[X.index, 'vendor'].values
mg = vendor == 'GE'
Xg = X.loc[mg]; yg = y[mg]
log('GE subgroup n=%d  MGMT+ %d / - %d  features=%d' % (len(yg), yg.sum(), (1 - yg).sum(), Xg.shape[1]))


def inner_select(Xtr, ytr, max_feat=200, corr_th=0.90, max_keep=40):
    auc = {}
    for c in Xtr.columns:
        v = Xtr[c].values
        m = np.isfinite(v)
        if m.sum() < 10 or len(np.unique(v[m])) < 2: continue
        try: auc[c] = abs(roc_auc_score(ytr[m], v[m]) - 0.5)
        except Exception: pass
    order = sorted(auc, key=auc.get, reverse=True)[:max_feat]
    kept, rk = [], Xtr[order].rank()
    for c in order:
        if all(abs(rk[c].corr(rk[k])) <= corr_th for k in kept):
            kept.append(c)
        if len(kept) >= max_keep: break
    return kept


def cv_auc(X_, y_, seed):
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    p = np.zeros(len(y_))
    for tr, te in skf.split(X_, y_):
        cols = inner_select(X_.iloc[tr], y_[tr])
        if not cols:
            p[te] = 0.5; continue
        pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()),
                         ('m', LogisticRegression(max_iter=5000, solver='liblinear'))])
        gs = GridSearchCV(pipe, {'m__C': [0.01, 0.1, 1, 10]}, cv=3, scoring='roc_auc', n_jobs=1)
        gs.fit(X_.iloc[tr][cols], y_[tr])
        p[te] = gs.predict_proba(X_.iloc[te][cols])[:, 1]
    return float(roc_auc_score(y_, p))


# 实际 AUC (同管线, 同 SEED, logistic)
obs = cv_auc(Xg, yg, SEED)
log('observed GE logistic AUC = %.4f' % obs)

rng = np.random.default_rng(SEED)
null = []
for b in range(N_PERM):
    yp = rng.permutation(yg)
    a = cv_auc(Xg, yp, SEED + 1000 + b)
    null.append(a)
    log('  perm %3d AUC=%.4f' % (b + 1, a))
null = np.array(null)
res = {'n_perm': N_PERM, 'observed_auc': obs,
       'null_mean': float(null.mean()), 'null_sd': float(null.std()),
       'null_p2.5': float(np.percentile(null, 2.5)), 'null_p97.5': float(np.percentile(null, 97.5)),
       'null_min': float(null.min()), 'null_max': float(null.max()),
       'perm_p_two_sided': float((np.sum(np.abs(null - 0.5) >= abs(obs - 0.5)) + 1) / (N_PERM + 1)),
       'frac_null_below_observed': float((null <= obs).mean())}
log(json.dumps(res, indent=2))
json.dump(res, open(os.path.join(ANA, 'ge_perm_null.json'), 'w'), indent=2)
np.savetxt(os.path.join(ANA, 'ge_perm_null_values.csv'), null, header='perm_auc')
print('saved -> analysis/ge_perm_null.json')
