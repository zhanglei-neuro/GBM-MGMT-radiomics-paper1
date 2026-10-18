# -*- coding: utf-8 -*-
r"""
Paper 1 - 完整分析脚本
题目: Multivendor MRI radiomics for predicting MGMT promoter methylation in
      glioblastoma using fully automated deep-learning segmentation:
      a cautionary negative-result study

运行: D:\Anaconda\python.exe  (numpy 1.26.4 / sklearn 1.5.1 / xgboost 3.4.1 /
      shap 0.46.0 / neuroCombat 0.2.12)

产出: analysis/*.json, analysis/*.csv, figures/*.png
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import os, json, warnings, sys
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from sklearn.calibration import calibration_curve
import xgboost as xgb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
warnings.filterwarnings('ignore')

SRC   = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57'
OUT   = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA   = os.path.join(OUT, 'analysis')
FIG   = os.path.join(OUT, 'figures')
os.makedirs(ANA, exist_ok=True); os.makedirs(FIG, exist_ok=True)
SEED  = 42
REPEATS = 5          # 重复分层 5 折 CV 的重复次数
NSPLIT  = 5

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 120})

def log(*a):
    print(*a, flush=True)

# ============================ 1. 数据装载 ============================
log('='*70); log('[1] loading'); log('='*70)
df = pd.read_csv(os.path.join(SRC, 'scripts', 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
Xall = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
Xall.index = [to_id(p) for p in Xall.index]        # 数据键统一为 Pxxx
co = pd.read_csv(os.path.join(SRC, 'scripts', 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]
common = [p for p in co.index if p in Xall.index]
# 2026-10-09 用户确认: 排除 United Imaging 1.5T 单例 (P021, 见稿件 3.8)
#   理由: n=1 无法作平台层面推断; ComBat 单例 batch 无法估计批参数; 敏感性分析证实剔除后结论不变
# 按 vendor 判断而非姓名: 发布版数据的患者标识为 Pxxx, 按姓名排除会静默失效并改变 n
common = [p for p in common if str(co.loc[p, 'vendor']).strip() != 'UIH']
X = Xall.loc[common].copy(); y = co.loc[common, 'y'].values
vendor = co.loc[common, 'vendor'].values
age = pd.to_numeric(co.loc[common, 'age'], errors='coerce').values
male = co.loc[common, 'gender'].astype(str).str.strip().eq('M').astype(int).values
log('N=%d  MGMT+ %d / MGMT- %d' % (len(y), int(y.sum()), int((1-y).sum())))
log('vendors: %s' % dict(pd.Series(vendor).value_counts()))

# 标签无关全局预处理
miss = X.isna().mean(); X = X[miss[miss <= 0.20].index]
vv = X.var(); X = X[vv[vv > 1e-12].index]
log('features after missing/zero-variance filter: %d' % X.shape[1])
n_feat_all = X.shape[1]

# ============================ 2. 通用建模组件 ============================
def inner_select(Xtr, ytr, max_feat=200, corr_th=0.90, max_keep=40):
    """内层: 单变量AUC排序 -> 贪心去相关 -> 至多 max_keep 个特征"""
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
    rk = Xtr[order].rank()
    for c in order:
        ok = True
        for k in kept:
            if abs(rk[c].corr(rk[k])) > corr_th:
                ok = False; break
        if ok:
            kept.append(c)
        if len(kept) >= max_keep:
            break
    return kept

def make_model(name, seed=SEED):
    if name == 'logistic':
        return (LogisticRegression(max_iter=5000, solver='liblinear'),
                {'C': [0.01, 0.1, 1, 10]})
    if name == 'xgboost':
        return (xgb.XGBClassifier(n_estimators=100, max_depth=2, learning_rate=0.1,
                                  subsample=0.8, colsample_bytree=0.8,
                                  eval_metric='logloss', random_state=seed, n_jobs=2,
                                  tree_method='hist'),
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
    gs = GridSearchCV(pipe, {'m__' + k: v for k, v in grid.items()},
                      cv=3, scoring='roc_auc', n_jobs=1)
    gs.fit(Xtr[cols], ytr)
    return gs.predict_proba(Xte[cols])[:, 1], cols

def youden(y_, p_):
    fpr, tpr, thr = roc_curve(y_, p_)
    i = int(np.argmax(tpr - fpr))
    pred = (p_ >= thr[i]).astype(int)
    tp = ((pred == 1) & (y_ == 1)).sum(); fn = ((pred == 0) & (y_ == 1)).sum()
    tn = ((pred == 0) & (y_ == 0)).sum(); fp = ((pred == 1) & (y_ == 0)).sum()
    return {'threshold': float(thr[i]), 'sens': float(tp/max(tp+fn,1)),
            'spec': float(tn/max(tn+fp,1)), 'acc': float((tp+tn)/len(y_))}

def rep_nested(name, X_, y_, repeats=REPEATS, seed=SEED):
    """重复分层5折嵌套CV -> 返回各重复 pooled OOF 概率与 AUC 列表"""
    aucs, probs, feat_sets = [], [], []
    for r in range(repeats):
        skf = StratifiedKFold(n_splits=NSPLIT, shuffle=True, random_state=seed + r)
        p = np.zeros(len(y_))
        for tr, te in skf.split(X_, y_):
            pp, cols = fit_predict(name, X_.iloc[tr], y_[tr], X_.iloc[te], seed + r)
            p[te] = pp; feat_sets.append(cols)
        aucs.append(roc_auc_score(y_, p)); probs.append(p)
    return np.array(aucs), probs, feat_sets

def boot_ci(y_, p_, n=2000, seed=SEED):
    rng = np.random.default_rng(seed)
    a = []
    for _ in range(n):
        idx = rng.integers(0, len(y_), len(y_))
        if len(np.unique(y_[idx])) < 2:
            continue
        a.append(roc_auc_score(y_[idx], p_[idx]))
    return float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))

# ============================ 3. 主分析: 嵌套CV ============================
log('='*70); log('[3] main nested CV (no vendor pre-filtering), repeats=%d' % REPEATS); log('='*70)
_cache_m = os.path.join(ANA, 'metrics_main.json')
_cache_p = os.path.join(ANA, 'predictions.csv')
main = {}
oof = pd.DataFrame({'patient': X.index, 'y': y, 'vendor': vendor})
if os.environ.get('REUSE_MAIN', '0') == '1' and os.path.exists(_cache_m) and os.path.exists(_cache_p):
    main = json.load(open(_cache_m))
    _cc = pd.read_csv(_cache_p, encoding='utf-8-sig')
    _cc['patient'] = [to_id(p) for p in _cc['patient']]   # 落盘表统一归一化为 Pxxx
    oof = _cc.set_index('patient').loc[X.index]
    oof = oof.reset_index()
    log('  [reuse] main CV loaded from cache: %s' % list(main.keys()))
else:
    for name in ['logistic', 'xgboost', 'randomforest']:
        aucs, probs, fsets = rep_nested(name, X, y)
        pmean = np.mean(probs, axis=0)
        lo, hi = boot_ci(y, pmean)
        met = youden(y, pmean)
        met.update({'auc': float(np.mean(aucs)), 'auc_sd': float(np.std(aucs)),
                    'auc_ci': [lo, hi], 'brier': float(brier_score_loss(y, pmean)),
                    'n_feat_selected_median': int(np.median([len(s) for s in fsets]))})
        main[name] = met
        oof['prob_' + name] = pmean
        log('  %-13s AUC=%.3f +/- %.3f (95%%CI %.3f-%.3f)  Acc=%.3f' %
            (name, met['auc'], met['auc_sd'], lo, hi, met['acc']))

    # 临床基线
    Xcl = pd.DataFrame({'age': age, 'male': male}, index=X.index)
    acls = []
    for r in range(REPEATS):
        skf = StratifiedKFold(n_splits=NSPLIT, shuffle=True, random_state=SEED + r)
        p = np.zeros(len(y))
        for tr, te in skf.split(Xcl, y):
            pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()),
                             ('m', LogisticRegression(max_iter=2000))])
            pipe.fit(Xcl.iloc[tr], y[tr]); p[te] = pipe.predict_proba(Xcl.iloc[te])[:, 1]
        acls.append(roc_auc_score(y, p)); oof['prob_clinical'] = p
    cl_lo, cl_hi = boot_ci(y, oof['prob_clinical'].values)
    main['clinical'] = {'auc': float(np.mean(acls)), 'auc_sd': float(np.std(acls)),
                        'auc_ci': [cl_lo, cl_hi],
                        'brier': float(brier_score_loss(y, oof['prob_clinical']))}
    log('  %-13s AUC=%.3f +/- %.3f (95%%CI %.3f-%.3f)' %
        ('clinical', main['clinical']['auc'], main['clinical']['auc_sd'], cl_lo, cl_hi))

    # 临床 + 影像 联合 (logistic)
    Xcomb = pd.concat([X, Xcl], axis=1)
    aucs_c, probs_c, _ = rep_nested('logistic', Xcomb, y)
    pc = np.mean(probs_c, axis=0); c_lo, c_hi = boot_ci(y, pc)
    main['combined'] = {'auc': float(np.mean(aucs_c)), 'auc_sd': float(np.std(aucs_c)),
                        'auc_ci': [c_lo, c_hi], 'brier': float(brier_score_loss(y, pc))}
    oof['prob_combined'] = pc
    log('  %-13s AUC=%.3f +/- %.3f (95%%CI %.3f-%.3f)' % ('combined', np.mean(aucs_c), np.std(aucs_c), c_lo, c_hi))

    json.dump(main, open(_cache_m, 'w'), indent=2)
    oof.to_csv(_cache_p, index=False, encoding='utf-8-sig')

# ============================ 4. 厂商鲁棒性 ============================
log('='*70); log('[4] vendor robustness: LOVO + within-vendor'); log('='*70)
ven = {}
# LOVO
lovo = {}
for v in ['Siemens', 'GE', 'UIH']:
    te = vendor == v; tr = ~te
    if te.sum() == 0:
        continue
    rec = {'n_test': int(te.sum())}
    for name in ['logistic', 'xgboost']:
        p, _ = fit_predict(name, X.iloc[tr], y[tr], X.iloc[te])
        rec[name + '_auc'] = float(roc_auc_score(y[te], p)) if len(np.unique(y[te])) > 1 else None
        oof.loc[te, 'lovo_' + name] = p
    if te.sum() >= 10:
        rec.update({('auc_' + k): v2 for k, v2 in youden(y[te], oof.loc[te, 'lovo_xgboost'].values).items()
                    if k in ('sens', 'spec', 'acc')})
    lovo[v] = rec
    lg, xg = rec.get('logistic_auc'), rec.get('xgboost_auc')
    log('  LOVO -> %-8s n=%2d  logistic=%s  xgboost=%s' %
        (v, te.sum(), ('%.3f' % lg) if lg is not None else 'NA',
         ('%.3f' % xg) if xg is not None else 'NA'))
ven['lovo'] = lovo

# 单厂商内部
within = {}
for v in ['Siemens', 'GE']:
    m = vendor == v
    if m.sum() < 25 or len(np.unique(y[m])) < 2:
        continue
    rec = {'n': int(m.sum()), 'n_pos': int(y[m].sum())}
    for name in ['logistic', 'xgboost']:
        a, pr, _ = rep_nested(name, X.loc[m], y[m])
        rec[name + '_auc'] = float(np.mean(a)); rec[name + '_sd'] = float(np.std(a))
        if name == 'logistic':
            oof.loc[m, 'within_logistic'] = np.mean(pr, axis=0)
    within[v] = rec
    log('  within %-8s n=%2d(+%d)  logistic=%.3f+/-%.3f  xgboost=%.3f+/-%.3f' %
        (v, rec['n'], rec['n_pos'], rec['logistic_auc'], rec['logistic_sd'],
         rec['xgboost_auc'], rec['xgboost_sd']))
ven['within'] = within
json.dump(ven, open(os.path.join(ANA, 'metrics_vendor.json'), 'w'), indent=2)
oof.to_csv(os.path.join(ANA, 'predictions.csv'), index=False, encoding='utf-8-sig')

# ============================ 5. 特征层面厂商差异 ============================
log('='*70); log('[5] feature-level vendor differences (KW + BH-FDR)'); log('='*70)
ms, mg = vendor == 'Siemens', vendor == 'GE'
kw_p, kw_eff = {}, {}
for c in X.columns:
    a = X.loc[ms, c].dropna(); b = X.loc[mg, c].dropna()
    if len(a) > 5 and len(b) > 5:
        try:
            kw_p[c] = stats.kruskal(a, b).pvalue
            # 效应量 rank-biserial 近似: 1 - 2U/(n1*n2)
            u = stats.mannwhitneyu(a, b, alternative='two-sided').statistic
            kw_eff[c] = 1 - 2 * u / (len(a) * len(b))
        except Exception:
            pass
ps = pd.Series(kw_p).sort_values()
n = len(ps); qv = ps.values * n / (np.arange(n) + 1)
qv = np.minimum.accumulate(qv[::-1])[::-1]
qser = pd.Series(qv, index=ps.index)
fdf = pd.DataFrame({'feature': ps.index, 'p': ps.values, 'q_bh': qser.values,
                    'effect': [kw_eff[c] for c in ps.index]})
fdf['abs_effect'] = fdf['effect'].abs()
fdf.to_csv(os.path.join(ANA, 'feature_vendor_kw.csv'), index=False, encoding='utf-8-sig')
fdr = {'n_tested': int(n), 'raw_p05': int((ps < 0.05).sum()),
       'q05': int((qser < 0.05).sum()), 'q01': int((qser < 0.01).sum()),
       'median_abs_effect': float(fdf['abs_effect'].median()),
       'q05_pct': float(100 * (qser < 0.05).mean())}
log('  tested=%d  raw p<0.05=%d (%.1f%%)  BH q<0.05=%d (%.1f%%)  q<0.01=%d' %
    (n, fdr['raw_p05'], 100*fdr['raw_p05']/n, fdr['q05'], fdr['q05_pct'], fdr['q01']))
json.dump(fdr, open(os.path.join(ANA, 'metrics_feature_fdr.json'), 'w'), indent=2)

# ============================ 6. ROI x SEQ 分层 ============================
log('='*70); log('[6] ROI x SEQ stratified AUC'); log('='*70)
roi_res = {}
for roi in ['TumorCore', 'Necrosis', 'Edema', 'FLAIRChanges']:
    for seq in ['T1CE', 'FLAIR', 'T1w']:
        cols = [c for c in X.columns if c.startswith(roi + '|' + seq + '|')]
        if len(cols) < 5:
            continue
        a, _, _ = rep_nested('logistic', X[cols], y, repeats=3)
        roi_res[roi + '|' + seq] = {'n_feat': len(cols), 'auc': float(np.mean(a)),
                                    'sd': float(np.std(a))}
        log('  %-22s nfeat=%4d  AUC=%.3f+/-%.3f' % (roi+'|'+seq, len(cols), np.mean(a), np.std(a)))
json.dump(roi_res, open(os.path.join(ANA, 'metrics_roi.json'), 'w'), indent=2)

# ============================ 7. ComBat 调和 ============================
log('='*70); log('[7] ComBat harmonisation (neuroCombat)'); log('='*70)
combat_res = {}
try:
    from neuroCombat import neuroCombat, neuroCombatFromTraining
    keep = vendor != 'UIH'          # UIH n=1 无法估计批参数
    Xc = X.loc[keep].copy(); yc = y[keep]; vc = vendor[keep]

    def combat_pair(Xtr_, vtr_, Xte_, vte_):
        """训练折拟合 ComBat -> 用训练估计变换测试折(无参数泄漏)"""
        med = Xtr_.median()
        Atr = Xtr_.fillna(med)
        est = neuroCombat(dat=Atr.T.values.astype(np.float64),
                          covars=pd.DataFrame({'vendor': vtr_}), batch_col='vendor')
        Btr = pd.DataFrame(est['data'].T, columns=Xtr_.columns, index=Xtr_.index)
        Ate = Xte_.fillna(med)
        try:
            out = neuroCombatFromTraining(dat=Ate.T.values.astype(np.float64),
                                          batch=vte_, estimates=est['estimates'])
            Bte = pd.DataFrame(out['data'].T, columns=Xte_.columns, index=Xte_.index)
        except Exception:
            Bte = Ate
        return Btr, Bte

    # 7a. 嵌套CV (每折内拟合 -> 测试折用训练估计变换)
    def rep_nested_combat(name, X_, y_, v_, repeats=REPEATS):
        aucs, probs = [], []
        for r in range(repeats):
            skf = StratifiedKFold(n_splits=NSPLIT, shuffle=True, random_state=SEED + r)
            p = np.zeros(len(y_))
            for tr, te in skf.split(X_, y_):
                Btr, Bte = combat_pair(X_.iloc[tr], v_[tr], X_.iloc[te], v_[te])
                pp, _ = fit_predict(name, Btr, y_[tr], Bte)
                p[te] = pp
            aucs.append(roc_auc_score(y_, p)); probs.append(p)
        return np.array(aucs), probs
    for name in ['logistic', 'xgboost']:
        a, pr = rep_nested_combat(name, Xc, yc, vc)
        combat_res[name] = {'auc': float(np.mean(a)), 'sd': float(np.std(a))}
        log('  ComBat nested CV %-9s AUC=%.3f+/-%.3f' % (name, np.mean(a), np.std(a)))
    # 7b. 调和效果: 调和后特征层面厂商差异
    med_all = Xc.median()
    est_all = neuroCombat(dat=Xc.fillna(med_all).T.values.astype(np.float64),
                          covars=pd.DataFrame({'vendor': vc}), batch_col='vendor')
    Xh = pd.DataFrame(est_all['data'].T, columns=Xc.columns, index=Xc.index)
    kw_after = 0
    for c in Xc.columns:
        a = Xh.loc[vc == 'Siemens', c].dropna(); b = Xh.loc[vc == 'GE', c].dropna()
        if len(a) > 5 and len(b) > 5 and stats.kruskal(a, b).pvalue < 0.05:
            kw_after += 1
    combat_res['kw_after'] = {'significant': int(kw_after), 'total': int(Xc.shape[1]),
                              'before_raw': int((ps < 0.05).sum()),
                              'before_q05': int((qser < 0.05).sum())}
    log('  after ComBat: vendor-significant features %d / %d (before raw %d / FDR %d)' %
        (kw_after, Xc.shape[1], int((ps < 0.05).sum()), int((qser < 0.05).sum())))
    # 7c. ComBat 后 LOVO (用全局估计: ComBat 固有地需要目标批数据 -> 需说明)
    lovo_h = {}
    for v in ['Siemens', 'GE']:
        te = vc == v; tr = ~te
        Btr = pd.DataFrame(est_all['data'].T, columns=Xc.columns, index=Xc.index).iloc[tr]
        out = neuroCombatFromTraining(dat=Xc.iloc[te].fillna(med_all).T.values.astype(np.float64),
                                      batch=vc[te], estimates=est_all['estimates'])
        Bte = pd.DataFrame(out['data'].T, columns=Xc.columns, index=Xc.index[te])
        p, _ = fit_predict('xgboost', Btr, yc[tr], Bte)
        lovo_h[v] = {'n_test': int(te.sum()), 'auc': float(roc_auc_score(yc[te], p))}
        log('  ComBat LOVO -> %-8s n=%2d AUC=%.3f' % (v, te.sum(), lovo_h[v]['auc']))
    combat_res['lovo'] = lovo_h
    # 7d. per-vendor z-score (VWS) 对照
    vws = {}
    for v in ['Siemens', 'GE']:
        te = vc == v; tr = ~te
        mu = Xc.iloc[tr].mean(); sd_ = Xc.iloc[tr].std().replace(0, 1)
        Xtr_z = (Xc.iloc[tr].fillna(mu) - mu) / sd_
        Xte_z = (Xc.iloc[te].fillna(Xc.iloc[te].mean()) - Xc.iloc[te].mean()) / Xc.iloc[te].std().replace(0, 1)
        p, _ = fit_predict('xgboost', Xtr_z, yc[tr], Xte_z)
        vws[v] = {'n_test': int(te.sum()), 'auc': float(roc_auc_score(yc[te], p))}
        log('  VWS LOVO -> %-8s n=%2d AUC=%.3f' % (v, te.sum(), vws[v]['auc']))
    combat_res['vws_lovo'] = vws
except Exception as e:
    import traceback; traceback.print_exc()
    log('  ComBat FAILED: %s' % e)
    combat_res['error'] = str(e)
json.dump(combat_res, open(os.path.join(ANA, 'metrics_combat.json'), 'w'), indent=2)

log('ALL ANALYSIS DONE')
