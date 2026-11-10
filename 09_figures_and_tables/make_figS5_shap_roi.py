# -*- coding: utf-8 -*-
r"""Figure S5 - exploratory attribution and stratification (layout-fixed version)

Rebuilds the exploratory composite from the deposited Zenodo intermediates only:

    (a) SHAP importance, XGBoost fitted on the full dataset (exploratory,
        not cross-validated) - recomputed here from
        zenodo/01_feature_matrix/radiomics_features_long.csv + the labels in
        zenodo/02_predictions_and_volumes/predictions_oof_main.csv;
    (b) AUC by ROI x sequence - zenodo/03_intermediate_results/metrics_roi.json;
    (c) main ROI x sequence combinations - same JSON.

Layout fix (2026-11-09, rev 2): (a) spans the full width; panels (b) and (c)
share the bottom row equally, with (b)'s leftmost text aligned to the leftmost
text of (a)'s y-tick labels and (c)'s right edge on (a)'s right axes edge
(measured at draw time, not hard-coded).

Run: D:/Anaconda/python.exe make_figS5_shap_roi.py
Output: <package>/supplementary_figures/FigureS5_exploratory_SHAP_ROI_sequence.png
"""

import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..', 'tools'))
from patient_ids import to_id  # noqa: E402  idempotent on P-coded data

import os, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
from xgboost import XGBClassifier
import xgboost as xgb

warnings.filterwarnings('ignore')

ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
ZEN = os.path.join(ROOT, 'zenodo')
OUT = os.path.join(ROOT, 'supplementary_figures',
                   'FigureS5_exploratory_SHAP_ROI_sequence.png')

SEED = 42
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5,
                     'axes.titlesize': 8.5, 'axes.labelsize': 7.5,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 300, 'figure.dpi': 150,
                     'savefig.bbox': 'tight'})

# ---------- data ------------------------------------------------------------
df = pd.read_csv(os.path.join(ZEN, '01_feature_matrix',
                             'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
X.index = [to_id(p) for p in X.index]

lab = pd.read_csv(os.path.join(ZEN, '02_predictions_and_volumes',
                               'predictions_oof_main.csv'), encoding='utf-8-sig')
lab['patient'] = [to_id(p) for p in lab['patient']]
lab = lab.set_index('patient')
common = [p for p in lab.index if p in X.index]
X = X.loc[common]
y = lab.loc[common, 'y'].astype(int).values

X = X[X.columns[X.isna().mean() <= 0.20]]
X = X[X.columns[X.var() > 1e-12]]
print('n =', len(y), ' features =', X.shape[1])

ROI = json.load(open(os.path.join(ZEN, '03_intermediate_results',
                                  'metrics_roi.json')))


# ---------- (a) SHAP: identical model to the published figure ----------------
def inner_select(Xtr, ytr, max_feat=200, corr_th=.90, max_keep=40):
    auc = {}
    for c in Xtr.columns:
        v = Xtr[c].values
        m = np.isfinite(v)
        if m.sum() < 10 or len(np.unique(v[m])) < 2:
            continue
        try:
            auc[c] = abs(roc_auc_score(ytr[m], v[m]) - .5)
        except Exception:
            pass
    order = sorted(auc, key=auc.get, reverse=True)[:max_feat]
    kept = []
    rk = Xtr[order].rank()
    for c in order:
        if all(abs(rk[c].corr(rk[k])) <= corr_th for k in kept):
            kept.append(c)
        if len(kept) >= max_keep:
            break
    return kept


cols = inner_select(X, y)
pipe = Pipeline([('i', SimpleImputer(strategy='median')),
                 ('s', StandardScaler()),
                 ('m', XGBClassifier(n_estimators=100, max_depth=2,
                                     learning_rate=.1, subsample=.8,
                                     colsample_bytree=.8, eval_metric='logloss',
                                     random_state=SEED, n_jobs=2,
                                     tree_method='hist'))])
pipe.fit(X[cols], y)
Xs2 = pipe[:-1].transform(X[cols])
contrib = pipe[-1].get_booster().predict(xgb.DMatrix(Xs2), pred_contribs=True)
sv = contrib[:, :-1]
if sv.ndim == 3:
    sv = sv[1]
ma = np.abs(sv).mean(0)
top = sorted(zip(cols, ma), key=lambda t: -t[1])[:20][::-1]
print('top SHAP:', '%.3f' % top[-1][1], top[-1][0])


def _short(f):
    nm = f.split('|', 2)[2]
    return nm if len(nm) <= 40 else nm[:39] + '\u2026'


names = ['%s [%s|%s]' % (_short(f), f.split('|')[0][:4], f.split('|')[1])
         for f, _ in top]

# ---------- figure: explicit rects so (b)/(c) align with (a) -----------------
L, R = 0.315, 0.985          # axes left/right edges shared by all panels
CBW, PAD = 0.015, 0.02       # colour-bar width and its pads
TICK = 0.075                 # room for the colour-bar and (c) y-tick labels
GAPW = PAD + CBW + PAD + TICK   # total space between the two bottom panels
W = (R - L - GAPW) / 2       # equal bottom-panel width

fig = plt.figure(figsize=(7.0, 6.4))
axa = fig.add_axes([L, 0.585, R - L, 0.360])
axb = fig.add_axes([L, 0.100, W, 0.315])
cax = fig.add_axes([L + W + PAD, 0.100, CBW, 0.315])
axc = fig.add_axes([L + W + GAPW, 0.100, W, 0.315])

# (a) SHAP
axa.barh(range(len(top)), [t[1] for t in top], color='#4C72B0')
axa.set_yticks(range(len(top)))
axa.set_yticklabels(names, fontsize=5.8)
axa.set_xlabel('mean |SHAP value|')
axa.set_title('(a) SHAP importance (XGBoost, full-data fit; exploratory, not cross-validated)',
              fontweight='bold', loc='left', pad=6)

# (b) heat map
rois = ['TumorCore', 'Necrosis', 'Edema', 'FLAIRChanges']
seqs = ['T1CE', 'FLAIR', 'T1w']
mat = np.full((len(rois), len(seqs)), np.nan)
for i, r in enumerate(rois):
    for j, s in enumerate(seqs):
        k = r + '|' + s
        if k in ROI:
            mat[i, j] = ROI[k]['auc']
im = axb.imshow(mat, cmap='RdBu_r', vmin=.4, vmax=.7, aspect='auto')
# aspect='auto' is required: with the default 'equal', apply_aspect() shrinks
# the axes box back to the 4x3 data aspect at every draw and breaks alignment.
for i in range(len(rois)):
    for j in range(len(seqs)):
        if not np.isnan(mat[i, j]):
            axb.text(j, i, '%.3f' % mat[i, j], ha='center', va='center',
                     fontsize=7, color='k', fontweight='bold')
axb.set_xticks(range(len(seqs)))
axb.set_xticklabels(seqs)
axb.set_yticks(range(len(rois)))
axb.set_yticklabels(rois, fontsize=6.5)
axb.set_title('(b) AUC by ROI \u00d7 sequence', fontweight='bold', loc='left', pad=6)
cb = fig.colorbar(im, cax=cax)
cb.set_ticks([.40, .50, .60, .70])
cb.ax.tick_params(labelsize=6)

# (c) main combinations
keys = ['TumorCore|T1CE', 'Necrosis|T1CE', 'Edema|FLAIR', 'FLAIRChanges|FLAIR']
keys = [k for k in keys if k in ROI]
v = [ROI[k]['auc'] for k in keys]
e = [ROI[k]['sd'] for k in keys]
bars = axc.bar(range(len(keys)), v, yerr=e, capsize=3, width=.62,
               color=['#4C72B0', '#DD8452', '#55A868', '#C44E52'])
for bb, a, ee in zip(bars, v, e):
    axc.text(bb.get_x() + bb.get_width() / 2, a + (ee or 0) + .012,
             '%.3f' % a, ha='center', fontsize=6.5)
axc.axhline(.5, color='k', ls='--', lw=.7)
axc.set_xticks(range(len(keys)))
axc.set_xticklabels([k.replace('|', '\n') for k in keys], fontsize=6.5)
axc.set_ylim(.4, .75)
axc.set_ylabel('AUC')
axc.set_title('(c) Main ROI \u00d7 sequence combinations',
              fontweight='bold', loc='left', pad=6)

os.makedirs(os.path.dirname(OUT), exist_ok=True)

# ---------- two-pass alignment of the bottom row -----------------------------
# Pass 1: draw once, measure where (a)'s y-tick text actually starts; move (b)
# so its own leftmost text lands exactly there, then give (b)/(c) equal widths
# spanning from that text edge to (a)'s right axes edge.
fig.canvas.draw()
rend = fig.canvas.get_renderer()
Wf = fig.bbox.width
a_left = min(t.get_window_extent(rend).x0 for t in axa.get_yticklabels()) / Wf
b_x0 = min(t.get_window_extent(rend).x0 for t in axb.get_yticklabels()) / Wf
B = L + (a_left - b_x0)                 # new left edge of (b)'s axes
y0, hgt = axb.get_position().y0, axb.get_position().height
W = (R - B - GAPW) / 2                  # equal widths; (c) right edge stays at R
axb.set_position([B, y0, W, hgt])
cax.set_position([B + W + PAD, y0, CBW, hgt])
axc.set_position([B + W + GAPW, y0, W, hgt])
print('aligned: (a) text left=%.3f -> (b) axes at %.3f, bottom width each %.3f'
      % (a_left, B, W))

fig.savefig(OUT)
print('saved', OUT)
