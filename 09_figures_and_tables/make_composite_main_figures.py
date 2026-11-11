# -*- coding: utf-8 -*-
r"""Paper 1 - composite main figures (Figure 3 of the six-figure set, and ESM Figure S5)

The submission uses six main figures. Four of them are the single-purpose panels produced by
`make_main_figures.py` (Figures 1, 3, 4 and 6 of the new numbering). The two composites built here
combine panels that are read together:

    Figure 3  (a) ROC curves   (b) calibration   (c) decision curve analysis
    Fig. S5   (a) SHAP importance   (b) AUC heat map by ROI x sequence   (c) main combinations
              (exploratory; moved out of the main set for the European Radiology version)

Panels are redrawn from the same analysis products as `make_main_figures.py`, at print width
(7 in) so that the 7-8.5 pt type is not scaled down when the figure is placed on the page.

Run: D:\Anaconda\python.exe make_composite_main_figures.py
Output: paper1/figures/Figure3_discrimination_ROC_cal_DCA.png
        paper1/supplementary_figures/FigureS5_exploratory_SHAP_ROI_sequence.png
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  patient code <-> local name
import os, json, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression  # noqa: F401  (kept for parity with make_main_figures)
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.calibration import calibration_curve
import xgboost as xgb
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis'); FIG = os.path.join(OUT, 'figures')
SUPP = os.path.join(OUT, 'supplementary_figures')
os.makedirs(FIG, exist_ok=True)
SEED = 42
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5, 'axes.titlesize': 8.5,
                     'axes.labelsize': 7.5, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 6.5, 'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 150, 'savefig.bbox': 'tight'})
CV = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3']

def log(*a): print(*a, flush=True)

# ---------- data ----------
df = pd.read_csv(os.path.join(SRC, 'scripts', 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
co = pd.read_csv(os.path.join(SRC, 'scripts', 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]
X.index = [to_id(p) for p in X.index]
common = [p for p in co.index if p in X.index]
common = [p for p in common if str(co.loc[p, 'vendor']).strip() != 'UIH']
X = X.loc[common]; y = co.loc[common, 'y'].values
X = X[X.columns[X.isna().mean() <= 0.20]]; X = X[X.columns[X.var() > 1e-12]]
pred = pd.read_csv(os.path.join(ANA, 'predictions.csv'), encoding='utf-8-sig')
pred['patient'] = [to_id(p) for p in pred['patient']]
pred = pred.set_index('patient').loc[X.index].reset_index()
ROI = json.load(open(os.path.join(ANA, 'metrics_roi.json')))
log('n =', len(y), 'features =', X.shape[1])


def dca(y_, p_, ths):
    n = len(y_); out = []
    for t in ths:
        pr = (p_ >= t).astype(int)
        tp = ((pr == 1) & (y_ == 1)).sum(); fp = ((pr == 1) & (y_ == 0)).sum()
        out.append(tp / n - fp / n * (t / (1 - t)))
    return out


# ============ Figure 2: discrimination, calibration, net benefit ============
log('Figure 2: ROC + calibration + DCA')
fig = plt.figure(figsize=(7.0, 5.9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.86], hspace=0.34, wspace=0.34,
                      width_ratios=[1.12, 0.88])
axa, axb, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

# (a) ROC
roc_lines = []
for name, lab, c in [('logistic', 'Radiomics — logistic', CV[0]),
                     ('xgboost', 'Radiomics — XGBoost', CV[1]),
                     ('randomforest', 'Radiomics — random forest', CV[2]),
                     ('clinical', 'Clinical (age + sex)', CV[3]),
                     ('combined', 'Clinical + radiomics', CV[4])]:
    p = pred['prob_' + name].values
    fpr, tpr, _ = roc_curve(y, p); a = roc_auc_score(y, p)
    roc_lines.append(axa.plot(fpr, tpr, color=c, lw=1.3, label='%s  %.3f' % (lab, a))[0])
axa.plot([0, 1], [0, 1], 'k--', lw=.7)
axa.set_xlabel('1 − Specificity'); axa.set_ylabel('Sensitivity')
axa.set_title('(a) ROC curves, pooled out-of-fold', fontweight='bold', loc='left')


def _legend_hits(ax, leg, lines, n=40):
    """Number of ROC segments that pass through the legend box (data coordinates)."""
    ax.figure.canvas.draw()
    bb = leg.get_window_extent(ax.figure.canvas.get_renderer()).transformed(ax.transData.inverted())
    x0, x1, y0, y1 = bb.x0, bb.x1, bb.y0, bb.y1
    hits = 0
    for ln in lines:
        xd, yd = np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)
        for i in range(len(xd) - 1):
            xs = np.linspace(xd[i], xd[i + 1], n)
            ys = np.linspace(yd[i], yd[i + 1], n)
            hits += int(((xs > x0) & (xs < x1) & (ys > y0) & (ys < y1)).any())
    return hits


def _roc_legend(ax, lines, frame=False):
    """Compact AUC legend in the lower-right corner, kept clear of the curves.

    The five ROC curves dip below the identity line in this cohort, so the empty corner is
    smaller than it looks. The legend is shrunk stepwise until no curve enters its bounding
    box; if even the smallest size collides, an opaque white box is used instead.
    """
    leg = None
    for fs, sp, hl in [(6.2, 0.24, 1.3), (5.8, 0.22, 1.2), (5.4, 0.20, 1.1), (5.0, 0.18, 1.0)]:
        if leg is not None:
            leg.remove()
        leg = ax.legend(loc='lower right', frameon=frame, fontsize=fs, labelspacing=sp,
                        handlelength=hl, handletextpad=0.35, borderpad=0.15,
                        borderaxespad=0.32, title='AUC', title_fontsize=fs + 0.2)
        if _legend_hits(ax, leg, lines) == 0:
            log('Fig. 3a legend: %.1f pt, clear of the curves' % fs)
            return leg
    leg.remove()
    if frame:
        return leg
    return _roc_legend(ax, lines, frame=True)


_roc_legend(axa, roc_lines)

# (b) calibration
for name, lab, c in [('logistic', 'Logistic', CV[0]), ('xgboost', 'XGBoost', CV[1]),
                     ('clinical', 'Clinical', CV[3])]:
    p = pred['prob_' + name].values
    fr, mn = calibration_curve(y, p, n_bins=5, strategy='quantile')
    axb.plot(mn, fr, 'o-', color=c, label=lab, ms=3.5, lw=1.2)
axb.plot([0, 1], [0, 1], 'k--', lw=.7)
axb.set_xlabel('Predicted probability'); axb.set_ylabel('Observed frequency')
axb.set_title('(b) Calibration', fontweight='bold', loc='left')
axb.legend(frameon=False, fontsize=6.5)

# (c) DCA
ths = np.linspace(.05, .9, 60)
for name, lab, c in [('logistic', 'Logistic', CV[0]), ('xgboost', 'XGBoost', CV[1])]:
    axc.plot(ths, dca(y, pred['prob_' + name].values, ths), color=c, lw=1.3, label=lab)
prev = y.mean()
axc.plot(ths, [prev - (1 - prev) * (t / (1 - t)) for t in ths], color='#999999', ls='--', lw=1.1,
         label='Treat all')
axc.axhline(0, color='k', lw=.7, ls=':')
axc.set_xlabel('Threshold probability'); axc.set_ylabel('Net benefit')
axc.set_ylim(-0.15, 0.5)
axc.set_title('(c) Decision curve analysis', fontweight='bold', loc='left')
axc.legend(frameon=False, fontsize=6.5, loc='upper right')
fig.savefig(os.path.join(FIG, 'Figure3_discrimination_ROC_cal_DCA.png'))
plt.close()

# ============ Figure 5: exploratory attribution and stratification ============
log('Figure 5: SHAP + ROI x sequence')


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
contrib = pipe[-1].get_booster().predict(xgb.DMatrix(Xs2), pred_contribs=True)
sv = contrib[:, :-1]
if sv.ndim == 3: sv = sv[1]
ma = np.abs(sv).mean(0)
top = sorted(zip(cols, ma), key=lambda t: -t[1])[:20][::-1]
def _short(f):
    nm = f.split('|', 2)[2]
    return nm if len(nm) <= 40 else nm[:39] + '\u2026'


names = ['%s [%s|%s]' % (_short(f), f.split('|')[0][:4], f.split('|')[1]) for f, _ in top]

fig = plt.figure(figsize=(7.0, 6.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], hspace=0.30, wspace=0.55)
axa, axb, axc = fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

# (a) SHAP
axa.barh(range(len(top)), [t[1] for t in top], color='#4C72B0')
axa.set_yticks(range(len(top))); axa.set_yticklabels(names, fontsize=5.8)
axa.set_xlabel('mean |SHAP value|')
axa.set_title('(a) SHAP importance (XGBoost, full-data fit; exploratory, not cross-validated)',
              fontweight='bold', loc='left')

# (b) heat map
rois = ['TumorCore', 'Necrosis', 'Edema', 'FLAIRChanges']; seqs = ['T1CE', 'FLAIR', 'T1w']
mat = np.full((len(rois), len(seqs)), np.nan)
for i, r in enumerate(rois):
    for j, s in enumerate(seqs):
        k = r + '|' + s
        if k in ROI: mat[i, j] = ROI[k]['auc']
im = axb.imshow(mat, cmap='RdBu_r', vmin=.4, vmax=.7)
for i in range(len(rois)):
    for j in range(len(seqs)):
        if not np.isnan(mat[i, j]):
            axb.text(j, i, '%.3f' % mat[i, j], ha='center', va='center', fontsize=7,
                     color='k', fontweight='bold')
axb.set_xticks(range(len(seqs))); axb.set_xticklabels(seqs)
axb.set_yticks(range(len(rois))); axb.set_yticklabels(rois, fontsize=6.5)
axb.set_title('(b) AUC by ROI × sequence', fontweight='bold', loc='left')
plt.colorbar(im, ax=axb, fraction=.05, pad=.02)

# (c) main combinations
keys = ['TumorCore|T1CE', 'Necrosis|T1CE', 'Edema|FLAIR', 'FLAIRChanges|FLAIR']
keys = [k for k in keys if k in ROI]
v = [ROI[k]['auc'] for k in keys]; e = [ROI[k]['sd'] for k in keys]
b = axc.bar(range(len(keys)), v, yerr=e, capsize=3, width=.62,
            color=['#4C72B0', '#DD8452', '#55A868', '#C44E52'])
for bb, a, ee in zip(b, v, e):
    axc.text(bb.get_x() + bb.get_width() / 2, a + (ee or 0) + .012, '%.3f' % a, ha='center', fontsize=6.5)
axc.axhline(.5, color='k', ls='--', lw=.7)
axc.set_xticks(range(len(keys)))
axc.set_xticklabels([k.replace('|', '\n') for k in keys], fontsize=6.5)
axc.set_ylim(.4, .75); axc.set_ylabel('AUC')
axc.set_title('(c) Main ROI × sequence combinations', fontweight='bold', loc='left')

fig.savefig(os.path.join(SUPP, 'FigureS5_exploratory_SHAP_ROI_sequence.png'))
plt.close()

log('COMPOSITE FIGURES DONE -> %s' % FIG)
