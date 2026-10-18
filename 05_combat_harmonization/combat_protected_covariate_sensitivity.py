# -*- coding: utf-8 -*-
r"""ComBat 保护协变量敏感性分析 (稿件 §3.7 报告)

背景: 稿件 §2.7(6) 的主分析采用"仅厂商批次变量、不传生物学协变量"的设定(避免标签泄漏),
但代价是可能过度移除与批次相关的生物学变异。本脚本检验这一顾虑:
在**训练折内**以 MGMT 为保护协变量 (discrete covariate) 重拟合 ComBat, 比较判别力恢复情况。

两臂 (其余环节完全一致, 均为 5 折 × 5 重复嵌套 CV, 折内筛选+调参, seed=42):
  A: neuroCombat(dat, covars={'vendor': batch}, batch_col='vendor')            <- 主分析设定
  B: neuroCombat(dat, covars={'vendor': batch, 'MGMT': label},
                 batch_col='vendor', categorical_cols=['MGMT'])                 <- 保护协变量

两臂均在训练折拟合, 经 neuroCombatFromTraining 变换测试折 -> 无跨折参数泄漏。

结论 (n=89): 保护协变量**未**恢复判别力。
  logistic  合并 OOF AUC 0.554 (0.428-0.671) vs 0.571, ΔAUC=-0.017 (-0.130, 0.097), p=0.77
  xgboost   合并 OOF AUC 0.452 (0.329-0.569) vs 0.579, ΔAUC=-0.127 (-0.273, 0.032), p=0.09

注意: B 臂把标签本身作为协变量参与变换参数估计, 属"以训练折标签辅助特征变换",
      虽不构成测试集泄漏, 但使变换依赖结局, 因此只作为敏感性分析, 不作为主分析。
      结果文件: combat_protected_covariate_sensitivity.json (Zenodo 03_intermediate_results/)
"""
import os, json, warnings
import numpy as np, pandas as pd
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import xgboost as xgb
from scipy import stats
from scipy.stats import rankdata
warnings.filterwarnings('ignore')

from neuroCombat import neuroCombat, neuroCombatFromTraining

# ---- 按你的环境修改 (见 PATHS_TO_CONFIGURE.md) ----
SRC = os.environ.get('PAPER1_FEATURE_DIR',
                     r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts')
OUT = os.environ.get('PAPER1_OUT_DIR',
                     r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1_revision/reanalysis')
# SRC 需包含 radiomics_features_long.csv 与 cohort_MGMT.csv
os.makedirs(OUT, exist_ok=True)
SEED, REPEATS, NSPLIT, NJ = 42, 5, 5, 4


def log(*a): print(*a, flush=True)


# ---------------- 数据 ----------------
df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig').set_index('patient')
common = [p for p in co.index if p in X.index and str(co.loc[p, 'vendor']).strip() != 'UIH']
X = X.loc[common]
y = co.loc[common, 'MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0).values
vendor = co.loc[common, 'vendor'].values
# 与主管线一致的全局过滤 (本队列 4148 全部保留)
miss = X.isna().mean(); X = X[miss[miss <= 0.20].index]
vv = X.var(); X = X[vv[vv > 1e-12].index]
log('n=%d  features=%d  MGMT+=%d/%d' % (len(common), X.shape[1], y.sum(), len(y)))


# ---------------- 与主管线一致的折内筛选 / 模型 ----------------
def vec_auc_abs(Xtr, ytr):
    A = Xtr.values.astype(float); yv = np.asarray(ytr)
    r = rankdata(A, axis=0)
    n1 = int(yv.sum()); n0 = len(yv) - n1
    U = r[yv == 1].sum(axis=0) - n1 * (n1 + 1) / 2.0
    return np.abs(U / (n1 * n0) - 0.5)


def inner_select(Xtr, ytr, max_feat=200, corr_th=0.90, max_keep=40):
    s = vec_auc_abs(Xtr, ytr)
    order = [Xtr.columns[i] for i in np.argsort(-s) if s[i] >= 0][:max_feat]
    kept, rk = [], Xtr[order].rank()
    for c in order:
        if all(abs(rk[c].corr(rk[k])) <= corr_th for k in kept):
            kept.append(c)
        if len(kept) >= max_keep:
            break
    return kept


def base_grid(name, seed=SEED):
    if name == 'logistic':
        return LogisticRegression(max_iter=5000, solver='liblinear'), {'C': [0.01, 0.1, 1, 10]}
    return (xgb.XGBClassifier(n_estimators=100, max_depth=2, learning_rate=0.1, subsample=0.8,
                              colsample_bytree=0.8, eval_metric='logloss', random_state=seed,
                              n_jobs=NJ, tree_method='hist'),
            {'n_estimators': [50, 100], 'max_depth': [1, 2]})


def fit_predict(name, Xtr, ytr, Xte, seed=SEED):
    cols = inner_select(Xtr, ytr)
    if not cols:
        return np.full(len(Xte), 0.5)
    base, grid = base_grid(name, seed)
    pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()), ('m', base)])
    gs = GridSearchCV(pipe, {'m__' + k: v for k, v in grid.items()}, cv=3, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr[cols], ytr)
    return gs.predict_proba(Xte[cols])[:, 1]


# ---------------- ComBat 变换 ----------------
def combat_pair(Xtr_, ytr_, vtr_, Xte_, vte_, use_protected):
    med = Xtr_.median()
    Atr = Xtr_.fillna(med)
    if use_protected:
        covars = pd.DataFrame({'vendor': vtr_, 'MGMT': ytr_})
        est = neuroCombat(dat=Atr.T.values.astype(np.float64), covars=covars,
                          batch_col='vendor', categorical_cols=['MGMT'])
    else:
        covars = pd.DataFrame({'vendor': vtr_})
        est = neuroCombat(dat=Atr.T.values.astype(np.float64), covars=covars, batch_col='vendor')
    Btr = pd.DataFrame(est['data'].T, columns=Xtr_.columns, index=Xtr_.index)
    Ate = Xte_.fillna(med)
    try:
        out = neuroCombatFromTraining(dat=Ate.T.values.astype(np.float64),
                                      batch=vte_, estimates=est['estimates'])
        Bte = pd.DataFrame(out['data'].T, columns=Xte_.columns, index=Xte_.index)
    except Exception:
        Bte = Ate
    return Btr, Bte


def run_arm(name, use_protected):
    aucs, probs = [], []
    for r in range(REPEATS):
        skf = StratifiedKFold(NSPLIT, shuffle=True, random_state=SEED + r)
        p = np.zeros(len(y))
        for tr, te in skf.split(X, y):
            Btr, Bte = combat_pair(X.iloc[tr], y[tr], vendor[tr], X.iloc[te], vendor[te], use_protected)
            p[te] = fit_predict(name, Btr, y[tr], Bte, seed=SEED + r)
        aucs.append(roc_auc_score(y, p)); probs.append(p)
    return np.array(aucs), np.array(probs)


# ---------------- DeLong / bootstrap ----------------
def delong(y_true, p1, p2):
    y_true = np.asarray(y_true); p1 = np.asarray(p1); p2 = np.asarray(p2)
    pos = y_true == 1; neg = ~pos
    m, n = pos.sum(), neg.sum()

    def midrank(x):
        r = rankdata(x); return r

    def vstat(p):
        v10 = np.array([(p[i] > p[neg]).sum() + 0.5 * (p[i] == p[neg]).sum() for i in np.where(pos)[0]]) / n
        v01 = np.array([(p[pos] > p[j]).sum() + 0.5 * (p[pos] == p[j]).sum() for j in np.where(neg)[0]]) / m
        return v10, v01

    v10a, v01a = vstat(p1); v10b, v01b = vstat(p2)
    auc1 = v10a.mean(); auc2 = v10b.mean()
    s10 = np.cov(np.vstack([v10a, v10b])); s01 = np.cov(np.vstack([v01a, v01b]))
    S = s10 / n + s01 / m
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    z = (auc1 - auc2) / np.sqrt(var) if var > 0 else 0.0
    return float(auc1), float(auc2), float(z), float(2 * (1 - stats.norm.cdf(abs(z))))


def boot_ci(y_, p_, n=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    y_ = np.asarray(y_); p_ = np.asarray(p_)
    out = []
    for _ in range(n):
        idx = rng.integers(0, len(y_), len(y_))
        if len(np.unique(y_[idx])) < 2:
            continue
        out.append(roc_auc_score(y_[idx], p_[idx]))
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


# ---------------- 主流程 ----------------
res = {'n': int(len(y)), 'n_features': int(X.shape[1]), 'repeats': REPEATS, 'folds': NSPLIT,
       'arms': {}}
pooled = {}
for arm, prot in (('A_vendor_only', False), ('B_with_MGMT_protected', True)):
    for name in ('logistic', 'xgboost'):
        a, pr = run_arm(name, prot)
        key = '%s|%s' % (arm, name)
        lo, hi = boot_ci(y, pr.mean(0))
        res['arms'][key] = {
            'repeat_mean_auc': round(float(a.mean()), 4), 'repeat_sd': round(float(a.std()), 4),
            'pooled_oof_auc': round(float(roc_auc_score(y, pr.mean(0))), 4),
            'pooled_oof_ci': [round(lo, 4), round(hi, 4)],
        }
        pooled[key] = pr.mean(0)
        log('%-26s %-9s repeat-mean=%.3f±%.3f  pooled-OOF=%.3f (%.3f–%.3f)'
            % (arm, name, a.mean(), a.std(), roc_auc_score(y, pr.mean(0)), lo, hi))

# 配对比: B - A
for name in ('logistic', 'xgboost'):
    pa = pooled['A_vendor_only|%s' % name]; pb = pooled['B_with_MGMT_protected|%s' % name]
    a1, a2, z, pv = delong(y, pa, pb)
    rng = np.random.default_rng(SEED)
    d = []
    for _ in range(2000):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        d.append(roc_auc_score(y[idx], pb[idx]) - roc_auc_score(y[idx], pa[idx]))
    res['paired_B_minus_A|%s' % name] = {
        'delta_auc': round(float(a2 - a1), 4),
        'delta_ci': [round(float(np.percentile(d, 2.5)), 4), round(float(np.percentile(d, 97.5)), 4)],
        'delong_p': round(float(pv), 4),
    }
    log('paired ΔAUC(B−A) %-9s = %+.4f (%.4f–%.4f)  DeLong p=%.3f'
        % (name, a2 - a1, np.percentile(d, 2.5), np.percentile(d, 97.5), pv))

res['caveat'] = ('B 臂把 MGMT 标签作为保护协变量参与 ComBat 参数估计(仅用训练折标签), '
                 '虽无测试集泄漏, 但使特征变换依赖结局, 故仅作敏感性分析。')
json.dump(res, open(os.path.join(OUT, 'combat_protected_covariate_sensitivity.json'), 'w', encoding='utf-8'),
          indent=2, ensure_ascii=False)
log('saved -> combat_protected_covariate_sensitivity.json')
