# -*- coding: utf-8 -*-
r"""Paper 1 - 论文级图表生成 (读 analysis/ 产物)
运行: D:\Anaconda\python.exe paper1_plots.py
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import os, json, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.calibration import calibration_curve
import xgboost as xgb
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis'); FIG = os.path.join(OUT, 'figures')
os.makedirs(FIG, exist_ok=True)
SEED = 42
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 120, 'savefig.bbox': 'tight'})
CV = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3']
VC = {'Siemens': '#4C72B0', 'GE': '#DD8452', 'UIH': '#55A868'}

def log(*a): print(*a, flush=True)

# ---------- 数据 ----------
df = pd.read_csv(os.path.join(SRC, 'scripts', 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
co = pd.read_csv(os.path.join(SRC, 'scripts', 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]            # 数据键统一为 Pxxx
X.index = [to_id(p) for p in X.index]
common = [p for p in co.index if p in X.index]
# 排除 United Imaging 1.5T 单例 (见稿件 3.8)。按 vendor 判断而非姓名
# (发布版数据的患者标识为 Pxxx, 按姓名排除会静默失效并改变 n)
common = [p for p in common if str(co.loc[p, 'vendor']).strip() != 'UIH']
X = X.loc[common]; y = co.loc[common, 'y'].values; vendor = co.loc[common, 'vendor'].values
X = X[X.columns[X.isna().mean() <= 0.20]]; X = X[X.columns[X.var() > 1e-12]]
pred = pd.read_csv(os.path.join(ANA, 'predictions.csv'), encoding='utf-8-sig')
pred['patient'] = [to_id(p) for p in pred['patient']]        # 落盘表统一归一化为 Pxxx
pred = pred.set_index('patient').loc[X.index].reset_index()  # 与 X / y 显式同行序
M = json.load(open(os.path.join(ANA, 'metrics_main.json')))
VEN = json.load(open(os.path.join(ANA, 'metrics_vendor.json')))
ROI = json.load(open(os.path.join(ANA, 'metrics_roi.json')))
CB = json.load(open(os.path.join(ANA, 'metrics_combat.json'))) if os.path.exists(
    os.path.join(ANA, 'metrics_combat.json')) else {}
FDR = json.load(open(os.path.join(ANA, 'metrics_feature_fdr.json')))
fdf = pd.read_csv(os.path.join(ANA, 'feature_vendor_kw.csv'), encoding='utf-8-sig')

# ============ Fig 1: batch effect (PCA / t-SNE by vendor & label) ============
log('fig1 batch effect')
Xs = StandardScaler().fit_transform(SimpleImputer(strategy='median').fit_transform(X))
pcs = PCA(n_components=2, random_state=SEED).fit_transform(Xs)
ts = TSNE(n_components=2, perplexity=20, random_state=SEED, init='pca', learning_rate='auto').fit_transform(Xs)
fig, ax = plt.subplots(2, 2, figsize=(11, 9.5))
for j, (emb, ttl) in enumerate([(pcs, 'PCA'), (ts, 't-SNE')]):
    # 厂商着色
    for v in ['Siemens', 'GE', 'UIH']:
        m = vendor == v
        if m.sum() == 0: continue
        ax[0, j].scatter(emb[m, 0], emb[m, 1], s=42, c=VC[v], label='%s (n=%d)' % (v, m.sum()),
                         alpha=.8, edgecolor='w', linewidth=.5)
    ax[0, j].set_title('%s — coloured by scanner vendor' % ttl, fontweight='bold')
    ax[0, j].legend(frameon=False, fontsize=9)
    ax[0, j].set_xlabel('%s 1' % ttl); ax[0, j].set_ylabel('%s 2' % ttl)
    # 标签着色
    for lab, c, nm in [(0, '#999999', 'MGMT-unmethylated'), (1, '#C44E52', 'MGMT-methylated')]:
        m = y == lab
        ax[1, j].scatter(emb[m, 0], emb[m, 1], s=42, c=c, label='%s (n=%d)' % (nm, m.sum()),
                         alpha=.8, edgecolor='w', linewidth=.5)
    ax[1, j].set_title('%s — coloured by MGMT status' % ttl, fontweight='bold')
    ax[1, j].legend(frameon=False, fontsize=9)
    ax[1, j].set_xlabel('%s 1' % ttl); ax[1, j].set_ylabel('%s 2' % ttl)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig1_batch_effect.png')); plt.close()

# ============ Fig 2: ROC ============
log('fig2 ROC')
plt.figure(figsize=(6.4, 6))
for name, lab, c in [('logistic', 'Radiomics — logistic', CV[0]),
                     ('xgboost', 'Radiomics — XGBoost', CV[1]),
                     ('randomforest', 'Radiomics — random forest', CV[2]),
                     ('clinical', 'Clinical (age + sex)', CV[3]),
                     ('combined', 'Clinical + radiomics', CV[4])]:
    p = pred['prob_' + name].values
    fpr, tpr, _ = roc_curve(y, p); a = roc_auc_score(y, p)
    plt.plot(fpr, tpr, color=c, lw=1.8, label='%s  AUC=%.3f' % (lab, a))
plt.plot([0, 1], [0, 1], 'k--', lw=.8)
plt.xlabel('1 − Specificity'); plt.ylabel('Sensitivity')
plt.title('MGMT prediction — repeated stratified 5-fold nested CV', fontweight='bold')
plt.legend(loc='lower right', frameon=False, fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig2_roc.png')); plt.close()

# ============ Fig 3: calibration + DCA ============
log('fig3 calibration + DCA')
fig, ax = plt.subplots(1, 2, figsize=(11.5, 5))
for name, lab, c in [('logistic', 'Logistic', CV[0]), ('xgboost', 'XGBoost', CV[1]),
                     ('clinical', 'Clinical', CV[3])]:
    p = pred['prob_' + name].values
    fr, mn = calibration_curve(y, p, n_bins=5, strategy='quantile')
    ax[0].plot(mn, fr, 'o-', color=c, label=lab)
ax[0].plot([0, 1], [0, 1], 'k--', lw=.8)
ax[0].set_xlabel('Predicted probability'); ax[0].set_ylabel('Observed frequency')
ax[0].set_title('Calibration', fontweight='bold'); ax[0].legend(frameon=False)

def dca(y_, p_, ths):
    n = len(y_); out = []
    for t in ths:
        pr = (p_ >= t).astype(int)
        tp = ((pr == 1) & (y_ == 1)).sum(); fp = ((pr == 1) & (y_ == 0)).sum()
        out.append(tp / n - fp / n * (t / (1 - t)))
    return out
ths = np.linspace(.05, .9, 60)
for name, lab, c in [('logistic', 'Logistic', CV[0]), ('xgboost', 'XGBoost', CV[1])]:
    ax[1].plot(ths, dca(y, pred['prob_' + name].values, ths), color=c, label=lab)
prev = y.mean()
ax[1].plot(ths, [prev - (1 - prev) * (t / (1 - t)) for t in ths], color='#999999', ls='--', label='Treat all')
ax[1].axhline(0, color='k', lw=.8, ls=':')
ax[1].set_xlabel('Threshold probability'); ax[1].set_ylabel('Net benefit')
ax[1].set_ylim(-0.15, 0.5)
ax[1].set_title('Decision curve analysis', fontweight='bold'); ax[1].legend(frameon=False)
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig3_calibration_dca.png')); plt.close()

# ============ Fig 4: vendor robustness ============
log('fig4 vendor robustness')
fig, ax = plt.subplots(1, 3, figsize=(14, 4.6))
# (a) LOVO
vs = ['Siemens', 'GE']; aucs = [VEN['lovo'][v].get('xgboost_auc') for v in vs]
ns = [VEN['lovo'][v]['n_test'] for v in vs]
bars = ax[0].bar(vs, aucs, color=[VC[v] for v in vs], width=.6)
for b, a, n in zip(bars, aucs, ns):
    if a is not None: ax[0].text(b.get_x() + b.get_width()/2, a + .02, 'AUC=%.3f\nn=%d' % (a, n), ha='center', fontsize=9)
ax[0].axhline(.5, color='k', ls='--', lw=.8)
ax[0].set_ylim(0, .9); ax[0].set_ylabel('AUC'); ax[0].set_title('(a) Leave-one-vendor-out (XGBoost)', fontweight='bold')
# (b) within-vendor
vs2 = ['Siemens', 'GE']; w1 = [VEN['within'][v]['logistic_auc'] for v in vs2]; w2 = [VEN['within'][v]['xgboost_auc'] for v in vs2]
xp = np.arange(len(vs2)); bw = .35
ax[1].bar(xp - bw/2, w1, bw, label='Logistic', color=CV[0])
ax[1].bar(xp + bw/2, w2, bw, label='XGBoost', color=CV[1])
for i, (a, b) in enumerate(zip(w1, w2)):
    ax[1].text(i - bw/2, a + .02, '%.3f' % a, ha='center', fontsize=8)
    ax[1].text(i + bw/2, b + .02, '%.3f' % b, ha='center', fontsize=8)
ax[1].axhline(.5, color='k', ls='--', lw=.8)
ax[1].set_xticks(xp); ax[1].set_xticklabels(['%s (n=%d)' % (v, VEN['within'][v]['n']) for v in vs2])
ax[1].set_ylim(0, .9); ax[1].set_ylabel('AUC'); ax[1].set_title('(b) Within-vendor nested CV', fontweight='bold'); ax[1].legend(frameon=False, fontsize=9)
# (c) ComBat / VWS
labels, vals, cols = [], [], []
if 'lovo' in CB:
    for v in ['Siemens', 'GE']:
        labels.append('ComBat\nLOVO %s' % v); vals.append(CB['lovo'][v]['auc']); cols.append(VC[v])
if 'vws_lovo' in CB:
    for v in ['Siemens', 'GE']:
        labels.append('VWS\nLOVO %s' % v); vals.append(CB['vws_lovo'][v]['auc']); cols.append(VC[v])
b3 = ax[2].bar(range(len(vals)), vals, color=cols, width=.6)
for b, a in zip(b3, vals):
    ax[2].text(b.get_x() + b.get_width()/2, a + .02, '%.3f' % a, ha='center', fontsize=9)
ax[2].axhline(.5, color='k', ls='--', lw=.8)
ax[2].set_xticks(range(len(labels))); ax[2].set_xticklabels(labels, fontsize=8)
ax[2].set_ylim(0, .9); ax[2].set_ylabel('AUC'); ax[2].set_title('(c) After harmonisation', fontweight='bold')
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig4_vendor_robustness.png')); plt.close()

# ============ Fig 5: feature-level vendor differences ============
log('fig5 feature vendor differences')
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5))
ax[0].hist(-np.log10(fdf['p'].values + 1e-300), bins=50, color='#4C72B0', alpha=.85)
ax[0].axvline(-np.log10(0.05), color='#C44E52', ls='--', label='p = 0.05')
ax[0].set_xlabel('−log10(p)  Kruskal–Wallis (Siemens vs GE)')
ax[0].set_ylabel('Number of features')
ax[0].set_title('(a) Scanner-vendor effect across %d features\n%.1f%% significant after BH-FDR' %
                (FDR['n_tested'], FDR['q05_pct']), fontweight='bold')
ax[0].legend(frameon=False)
topf = fdf.nlargest(15, 'abs_effect').iloc[::-1]
short = [f.split('|', 2)[2][:34] + '  [' + f.split('|')[0][:4] + ']' for f in topf['feature']]
ax[1].barh(range(len(topf)), topf['abs_effect'].values, color='#DD8452')
ax[1].set_yticks(range(len(topf))); ax[1].set_yticklabels(short, fontsize=7.5)
ax[1].set_xlabel('|rank-biserial effect size|')
ax[1].set_title('(b) Top-15 features by vendor effect size', fontweight='bold')
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig5_feature_vendor.png')); plt.close()

# ============ Fig 6: SHAP ============
log('fig6 SHAP')
try:
    def inner_select(Xtr, ytr, max_feat=200, corr_th=.90, max_keep=40):
        auc = {}
        for c in Xtr.columns:
            v = Xtr[c].values; m = np.isfinite(v)
            if m.sum() < 10 or len(np.unique(v[m])) < 2: continue
            try: auc[c] = abs(roc_auc_score(ytr[m], v[m]) - .5)
            except Exception: pass
        order = sorted(auc, key=auc.get, reverse=True)[:max_feat]
        kept = []; rk = Xtr[order].rank()
        for c in order:
            if all(abs(rk[c].corr(rk[k])) <= corr_th for k in kept): kept.append(c)
            if len(kept) >= max_keep: break
        return kept
    cols = inner_select(X, y)
    pipe = Pipeline([('i', SimpleImputer(strategy='median')), ('s', StandardScaler()),
                     ('m', xgb.XGBClassifier(n_estimators=100, max_depth=2, learning_rate=.1,
                                             subsample=.8, colsample_bytree=.8, eval_metric='logloss',
                                             random_state=SEED, n_jobs=2, tree_method='hist'))])
    pipe.fit(X[cols], y)
    Xs2 = pipe[:-1].transform(X[cols])
    dm = xgb.DMatrix(Xs2)
    contrib = pipe[-1].get_booster().predict(dm, pred_contribs=True)
    sv = contrib[:, :-1]
    if sv.ndim == 3: sv = sv[1]
    ma = np.abs(sv).mean(0)
    top = sorted(zip(cols, ma), key=lambda t: -t[1])[:20][::-1]
    names = ['%s [%s|%s]' % (f.split('|', 2)[2][:28], f.split('|')[0][:4], f.split('|')[1]) for f, _ in top]
    plt.figure(figsize=(8.5, 7))
    plt.barh(range(len(top)), [t[1] for t in top], color='#4C72B0')
    plt.yticks(range(len(top)), names, fontsize=7.5)
    plt.xlabel('mean |SHAP value|')
    plt.title('SHAP feature importance (XGBoost, full-data fit)', fontweight='bold')
    plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig6_shap.png')); plt.close()
    pd.DataFrame(top, columns=['feature', 'mean_abs_shap']).to_csv(
        os.path.join(ANA, 'selected_features_shap.csv'), index=False, encoding='utf-8-sig')
except Exception as e:
    log('  SHAP failed:', e)

# ============ Fig 7: ROI x SEQ AUC ============
log('fig7 ROI x SEQ')
rois = ['TumorCore', 'Necrosis', 'Edema', 'FLAIRChanges']; seqs = ['T1CE', 'FLAIR', 'T1w']
mat = np.full((len(rois), len(seqs)), np.nan)
for i, r in enumerate(rois):
    for j, s in enumerate(seqs):
        k = r + '|' + s
        if k in ROI: mat[i, j] = ROI[k]['auc']
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.8))
im = ax[0].imshow(mat, cmap='RdBu_r', vmin=.4, vmax=.7)
for i in range(len(rois)):
    for j in range(len(seqs)):
        if not np.isnan(mat[i, j]):
            ax[0].text(j, i, '%.3f' % mat[i, j], ha='center', va='center', fontsize=10,
                       color='k', fontweight='bold')
ax[0].set_xticks(range(len(seqs))); ax[0].set_xticklabels(seqs)
ax[0].set_yticks(range(len(rois))); ax[0].set_yticklabels(rois)
ax[0].set_title('(a) AUC by ROI × sequence (nested CV)', fontweight='bold')
plt.colorbar(im, ax=ax[0], fraction=.046, label='AUC')
# 主 combo 柱状
keys = ['TumorCore|T1CE', 'Necrosis|T1CE', 'Edema|FLAIR', 'FLAIRChanges|FLAIR']
keys = [k for k in keys if k in ROI]
v = [ROI[k]['auc'] for k in keys]; e = [ROI[k]['sd'] for k in keys]
b = ax[1].bar(range(len(keys)), v, yerr=e, capsize=4, color=['#4C72B0', '#DD8452', '#55A868', '#C44E52'])
for bb, a in zip(b, v): ax[1].text(bb.get_x()+bb.get_width()/2, a+.015, '%.3f' % a, ha='center', fontsize=9)
ax[1].axhline(.5, color='k', ls='--', lw=.8)
ax[1].set_xticks(range(len(keys))); ax[1].set_xticklabels([k.replace('|', '\n') for k in keys], fontsize=8)
ax[1].set_ylim(.4, .75); ax[1].set_ylabel('AUC')
ax[1].set_title('(b) Primary ROI × sequence combinations', fontweight='bold')
plt.tight_layout(); plt.savefig(os.path.join(FIG, 'fig7_roi_auc.png')); plt.close()

log('ALL FIGURES DONE -> %s' % FIG)
