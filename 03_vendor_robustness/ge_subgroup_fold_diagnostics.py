# -*- coding: utf-8 -*-
r"""Paper 1 -- R1-L2: GE 亚组内 AUC<0.5 的可靠性诊断

回答审稿人: 报告 GE 亚组内 (a) 每折训练/测试的类别分布, (b) 每折 AUC 与方向,
(c) 内层筛选出的特征在训练折上的关联方向及其在测试折上的稳定性,
(d) 单变量层面 (不做特征选择) 的标签-特征关联方向统计。

输出: analysis/ge_diagnostics.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import os, json, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
ANA = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/analysis'
SEED = 42

df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
Xa = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
Xa.index = [to_id(p) for p in Xa.index]            # 数据键统一为 Pxxx
co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().str.startswith('+').astype(int)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]
# 排除 United Imaging 1.5T 单例 (见稿件 3.8)。按 vendor 判断, 不按姓名:
# 发布版数据的患者标识为 Pxxx, 按姓名排除会静默失效并改变 n。
common = [p for p in co.index if p in Xa.index
          and str(co.loc[p, 'vendor']).strip() != 'UIH']
X = Xa.loc[common]
X = X[X.columns[X.isna().mean() <= 0.20]]; X = X[X.columns[X.var() > 1e-12]]
y = co.loc[X.index, 'y'].values
mg = co.loc[X.index, 'vendor'].values == 'GE'
Xg, yg = X.loc[mg], y[mg]
n, npos = len(yg), int(yg.sum())
print('GE n=%d  MGMT+ %d / - %d' % (n, npos, n - npos))

def inner_select(Xtr, ytr, max_feat=200, corr_th=0.90, max_keep=40):
    auc = {}
    for c in Xtr.columns:
        v = Xtr[c].values; m = np.isfinite(v)
        if m.sum() < 10 or len(np.unique(v[m])) < 2: continue
        try: auc[c] = abs(roc_auc_score(ytr[m], v[m]) - 0.5)
        except Exception: pass
    order = sorted(auc, key=auc.get, reverse=True)[:max_feat]
    kept, rk = [], Xtr[order].rank()
    for c in order:
        if all(abs(rk[c].corr(rk[k])) <= corr_th for k in kept):
            kept.append(c)
        if len(kept) >= max_keep: break
    return kept, auc

res = {'n': n, 'n_pos': npos, 'n_neg': n - npos}

# ---- (a)(b) 每折类别分布与 AUC ----
skf = StratifiedKFold(5, shuffle=True, random_state=SEED)
folds = []
oof = np.zeros(n)
signs = []
for k, (tr, te) in enumerate(skf.split(Xg, yg)):
    cols, auc = inner_select(Xg.iloc[tr], yg[tr])
    pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()),
                     ('m', LogisticRegression(max_iter=5000, solver='liblinear'))])
    gs = GridSearchCV(pipe, {'m__C': [0.01, 0.1, 1, 10]}, cv=3, scoring='roc_auc', n_jobs=1)
    gs.fit(Xg.iloc[tr][cols], yg[tr])
    p = gs.predict_proba(Xg.iloc[te][cols])[:, 1]
    oof[te] = p
    try: fa = float(roc_auc_score(yg[te], p))
    except Exception: fa = np.nan
    # 训练折上 top10 特征的关联方向
    top = sorted(auc, key=auc.get, reverse=True)[:10]
    sd = []
    for c in top:
        v = Xg.iloc[tr][c].values; m = np.isfinite(v)
        if m.sum() < 10: continue
        try: sd.append(float(roc_auc_score(yg[tr][m], v[m])))       # >0.5 表示值高->阳性
        except Exception: pass
    folds.append({'fold': k + 1, 'n_train': int(len(tr)), 'n_test': int(len(te)),
                  'train_pos': int(yg[tr].sum()), 'train_neg': int((1 - yg[tr]).sum()),
                  'test_pos': int(yg[te].sum()), 'test_neg': int((1 - yg[te]).sum()),
                  'n_selected': len(cols), 'test_auc': fa})
    signs += sd
    print('fold %d  train %d (%d+/%d-)  test %d (%d+/%d-)  nsel=%d  testAUC=%.3f'
          % (k + 1, len(tr), yg[tr].sum(), (1 - yg[tr]).sum(), len(te), yg[te].sum(), (1 - yg[te]).sum(), len(cols), fa))
res['folds'] = folds
res['pooled_oof_auc'] = float(roc_auc_score(yg, oof))
print('pooled OOF AUC = %.4f' % res['pooled_oof_auc'])

sd = np.array(signs)
res['top10_direction_on_train'] = {
    'n_observed': int(len(sd)),
    'frac_auc_above_0.5': float((sd > 0.5).mean()),
    'median_auc': float(np.median(sd)),
    'pct_below_0.4': float((sd < 0.4).mean()),
    'pct_above_0.6': float((sd > 0.6).mean()),
    'note': '训练折内 top10 单变量 AUC (值高->阳性); 50% 左右说明方向随机'
}

# ---- (d) 全 GE 亚组单变量方向分布 ----
uni = []
for c in Xg.columns:
    v = Xg[c].values; m = np.isfinite(v)
    if m.sum() < 20 or len(np.unique(v[m])) < 2: continue
    try: uni.append(float(roc_auc_score(yg[m], v[m])))
    except Exception: pass
uni = np.array(uni)
res['univariate_all_features_GE'] = {
    'n': int(len(uni)), 'frac_above_0.5': float((uni > 0.5).mean()),
    'frac_signif_p05_by_perm_note': '见 ge_perm_null.json',
    'q10': float(np.percentile(uni, 10)), 'q50': float(np.percentile(uni, 50)),
    'q90': float(np.percentile(uni, 90)),
    'note': '若存在真实的系统性反向关联, 单变量 AUC 分布应整体左移 (frac_above_0.5 明显 <0.5)'
}
print('univariate frac>0.5 = %.3f  median = %.3f' % (res['univariate_all_features_GE']['frac_above_0.5'],
                                                    res['univariate_all_features_GE']['q50']))

json.dump(res, open(os.path.join(ANA, 'ge_diagnostics.json'), 'w'), indent=2)
print('saved -> analysis/ge_diagnostics.json')
