# -*- coding: utf-8 -*-
r"""Paper 1 -- 补充图 S3 (批次效应定量化) 与 S4 (厂商 vs 采集协议 去混杂分解)"""
import os, json, warnings
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

ANA = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/analysis'
FIG = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/figures'
os.makedirs(FIG, exist_ok=True)
plt.rcParams['font.size'] = 9
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# ---------------- 图 S3 ----------------
rev = json.load(open(os.path.join(ANA, 'review_metrics.json')))
asw = rev['asw']; pca = rev['pca_variance_between_groups']

fig, axes = plt.subplots(1, 3, figsize=(11, 3.2))
# (a) ASW 对比
ax = axes[0]
labels = ['Vendor\n(Siemens vs GE)', 'MGMT\n(methylated vs not)']
vals = [asw['vendor']['asw'], asw['mgmt']['asw']]
ps = [asw['vendor']['perm_p'], asw['mgmt']['perm_p']]
cols = ['#C00000', '#404040']
b = ax.bar(labels, vals, color=cols, width=0.55)
for r, v, p in zip(b, vals, ps):
    ax.text(r.get_x() + r.get_width() / 2, v + 0.006, '%.3f\n(p=%.3f)' % (v, p),
            ha='center', va='bottom', fontsize=8.5)
ax.set_ylabel('Average silhouette width (ASW)')
ax.set_title('(a) Quantitative batch separation', fontsize=9.5)
ax.set_ylim(0, max(vals) * 1.45)
ax.axhline(0, color='k', lw=0.6)

# (b) PC1-10 组间变异占比
ax = axes[1]
v2 = [pca['vendor_between_group_share_pc1_10'], pca['mgmt_between_group_share_pc1_10']]
b = ax.bar(['Vendor', 'MGMT'], [x * 100 for x in v2], color=cols, width=0.55)
for r, v in zip(b, v2):
    ax.text(r.get_x() + r.get_width() / 2, v * 100 + 0.8, '%.1f%%' % (v * 100),
            ha='center', va='bottom', fontsize=8.5)
ax.set_ylabel('Between-group variance share of PC1–10 (%)')
ax.set_title('(b) Variance decomposition', fontsize=9.5)
ax.set_ylim(0, max(v2) * 100 * 1.3)

# (c) PC 方差解释率 (scree)
ax = axes[2]
sc = pca['top10_pc_variance_explained']
ax.bar(range(1, len(sc) + 1), [x * 100 for x in sc], color='#4472C4', width=0.6)
ax.set_xlabel('Principal component'); ax.set_ylabel('Variance explained (%)')
ax.set_title('(c) PCA scree (top 10 PCs)', fontsize=9.5)
ax.set_xticks(range(1, len(sc) + 1))

fig.tight_layout()
fig.savefig(os.path.join(FIG, 'FigureS3_batch_effect_quantification.png'), dpi=300)
plt.close(fig)
print('saved FigureS3')

# ---------------- 图 S4 ----------------
m4 = json.load(open(os.path.join(ANA, 'm4_vendor_vs_protocol.json')))
fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4))

# (a) 显著特征数对比
ax = axes[0]
u = m4['unadjusted_vendor_only']['q05']
a = m4['adjusted_for_acquisition']['q05']
b = ax.bar(['Vendor only\n(Kruskal–Wallis)', 'Vendor adjusted for\nsequence type, thickness,\nTR / TE / TI'],
           [u, a], color=['#C00000', '#2E75B6'], width=0.5)
for r, v in zip(b, [u, a]):
    ax.text(r.get_x() + r.get_width() / 2, v + 60, '%d\n(%.1f%%)' % (v, 100 * v / 4148),
            ha='center', va='bottom', fontsize=9)
ax.set_ylabel('BH-FDR q < 0.05 features (out of 4148)')
ax.set_ylim(0, u * 1.28)
ax.set_title('(a) Vendor-significant features\nbefore vs after adjustment', fontsize=9.5)

# (b) 偏 R² 对比
ax = axes[1]
pv = m4['partial_r2']
vals = [pv['vendor_given_acquisition_median'], pv['acquisition_given_vendor_median']]
bb = ax.barh(['Vendor | acquisition\n(median partial R²)', 'Acquisition | vendor\n(median partial R²)'],
             vals, color=['#C00000', '#2E75B6'], height=0.5)
for r, v in zip(bb, vals):
    ax.text(v + 0.002, r.get_y() + r.get_height() / 2, '%.3f' % v, va='center', fontsize=9)
ax.set_xlabel('Partial R²'); ax.set_xlim(0, max(vals) * 1.4)
ax.set_title('(b) Which factor explains more feature variance?', fontsize=9.5)

fig.tight_layout()
fig.savefig(os.path.join(FIG, 'FigureS4_vendor_vs_protocol.png'), dpi=300)
plt.close(fig)
print('saved FigureS4')
