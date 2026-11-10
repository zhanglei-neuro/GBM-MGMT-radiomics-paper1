# -*- coding: utf-8 -*-
r"""Figure 6 - volumetric markers (x-axis fix for panel a)

Rebuilds the three panels of Figure 6 from the deposited Zenodo intermediates:

    (a) distribution of the total-tumour/oedema volume ratio (TT = TC + NE)
        by MGMT status - per_case_volumes.csv + predictions_oof_main.csv;
    (b) univariate AUCs of the seven markers - table3_main.csv;
    (c) TT/ED AUC by scanner vendor - recomputed from the per-case table.

Fix (2026-11-09): panel (a) previously hard-coded xlim(0, 6) while the largest
observed ratio is 1.47, which squeezed the histogram into a quarter of the
panel. The axis now spans the observed data range.

Run: D:/Anaconda/python.exe make_figure6_volume_biomarkers.py
Output: <package>/figures/Figure6_volume_biomarkers.png
"""

import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _os.path.join(_HERE, '..', 'tools'))

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

ROOT = os.path.normpath(os.path.join(_HERE, '..', '..'))
ZEN = os.path.join(ROOT, 'zenodo')
OUT = os.path.join(ROOT, 'figures', 'Figure6_volume_biomarkers.png')

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7.5,
                     'axes.titlesize': 8.5, 'axes.labelsize': 7.5,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 6.5, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 300,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})

# ---------- data ------------------------------------------------------------
v = pd.read_csv(os.path.join(ZEN, '02_predictions_and_volumes',
                             'per_case_volumes.csv'), encoding='utf-8-sig')
lab = pd.read_csv(os.path.join(ZEN, '02_predictions_and_volumes',
                               'predictions_oof_main.csv'), encoding='utf-8-sig')
d = v.merge(lab[['patient', 'y', 'vendor']], on='patient')
assert len(d) == 89, len(d)
d['ratio'] = (d.TC + d.NE) / d.ED

t3 = pd.read_csv(os.path.join(ZEN, '03_intermediate_results',
                              'table3_main.csv'), encoding='utf-8-sig')
auc = dict(zip(t3['Marker'], t3['AUC'].astype(float)))

ZH = 0.687  # Zahirovic et al. 2025 external benchmark

# ---------- figure ----------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.9))
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.20, top=0.86, wspace=0.42)

# (a) TT/ED ratio distribution - x-axis now spans the observed data
bins = np.linspace(0, 1.6, 19)
for val, nm, c in [(0, 'MGMT\u2212 (n=%d)' % (d.y == 0).sum(), '#999999'),
                   (1, 'MGMT+ (n=%d)' % (d.y == 1).sum(), '#C44E52')]:
    ax[0].hist(d.loc[d.y == val, 'ratio'], bins=bins, alpha=.6, color=c, label=nm)
ax[0].set_xlim(0, 1.6)
ax[0].set_xlabel('Total-tumour / oedema volume ratio')
ax[0].set_ylabel('Number of patients')
ax[0].legend(frameon=False, fontsize=6.5, loc='upper right',
             bbox_to_anchor=(0.99, 0.99), borderaxespad=0)
ax[0].set_title('(a) TT/ED ratio by MGMT status', fontweight='bold', loc='left', pad=6)

# (b) univariate AUCs of the seven markers
order = ['TumorCore', 'TotalTumor', 'TC_ED_ratio', 'TT_ED_ratio',
         'Edema', 'Necrosis', 'FLAIRChanges'][::-1]
vals = [auc[m] for m in order]
ax[1].barh(range(len(order)), vals, color='#4C72B0', height=.62)
for i, a_ in enumerate(vals):
    ax[1].text(a_ + .012, i, '%.3f' % a_, va='center', fontsize=6.5)
ax[1].set_yticks(range(len(order)))
ax[1].set_yticklabels(order, fontsize=6.5)
ax[1].axvline(.5, color='k', ls=':', lw=.7)
ax[1].axvline(ZH, color='#C44E52', ls=':', lw=.8)
ax[1].text(ZH + .01, 0.1, 'Zahirovic 2025 (%.3f)' % ZH, color='#C44E52',
           fontsize=6, ha='left', va='bottom')
ax[1].set_xlim(.3, .85)
ax[1].set_xlabel('Univariate AUC (higher value \u2192 MGMT+)')
ax[1].set_title('(b) Markers vs benchmark', fontweight='bold', loc='left', pad=6)

# (c) TT/ED AUC by vendor, recomputed from the per-case table
groups = [('All (n=89)', d), ('Siemens', d[d.vendor == 'Siemens']),
          ('GE', d[d.vendor == 'GE'])]
gvals = [roc_auc_score(g.y, g.ratio) for _, g in groups]
bars = ax[2].bar(range(3), gvals, width=.62,
                 color=['#8172B3', '#4C72B0', '#DD8452'])
for bb, a_ in zip(bars, gvals):
    ax[2].text(bb.get_x() + bb.get_width() / 2, a_ + .012, '%.3f' % a_,
               ha='center', fontsize=6.5)
ax[2].axhline(.5, color='k', ls='--', lw=.7)
ax[2].axhline(ZH, color='#C44E52', ls=':', lw=.8)
ax[2].text(2.45, .718, 'Zahirovic 2025 (%.3f)' % ZH, color='#C44E52',
           fontsize=6, ha='right', va='bottom')
ax[2].set_xticks(range(3))
ax[2].set_xticklabels([n for n, _ in groups])
ax[2].set_ylim(.3, .8)
ax[2].set_ylabel('AUC')
ax[2].set_title('(c) TT/ED AUC by vendor', fontweight='bold', loc='left', pad=6)
print('TT/ED AUC by vendor:', ['%.3f' % x for x in gvals])

os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT)
print('saved', OUT)
