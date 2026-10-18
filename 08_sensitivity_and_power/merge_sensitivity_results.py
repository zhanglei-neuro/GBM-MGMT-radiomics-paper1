# -*- coding: utf-8 -*-
r"""合并 4 个并行敏感性分析结果 -> analysis/sensitivity_cohorts.json + tables/S4_*.csv/.md

用法: python merge_sens.py
"""
import os
import json
import pandas as pd

OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
TB = os.path.join(OUT, 'supplementary', 'tables')

ORDER = ['A', 'B', 'C', 'D']
LABELS = {
    'A': 'A. Full cohort incl. the single UIH 1.5 T case',
    'B': 'B. PRIMARY cohort (UIH case excluded)',
    'C': 'C. 3D thin-slice T1c only',
    'D': 'D. 2D thick-slice T1c only',
}

res = []
for tag in ORDER:
    fp = os.path.join(ANA, 'sens_cohort_%s.json' % tag)
    r = json.load(open(fp, encoding='utf-8'))
    r['tag'] = tag
    res.append(r)

json.dump(res, open(os.path.join(ANA, 'sensitivity_cohorts.json'), 'w'),
          indent=2, ensure_ascii=False)

rows = []
for r in res:
    m = r['models']
    rows.append({
        'Cohort': LABELS[r['tag']],
        'n': r['n'],
        'MGMT+': r['n_pos'],
        'Logistic AUC': m['logistic']['auc'],
        'XGBoost AUC': m['xgboost']['auc'],
        'Random forest AUC': m['randomforest']['auc'],
        'Clinical AUC': m['clinical']['auc'],
        'Clinical 95% CI': '%.3f-%.3f' % tuple(m['clinical']['ci']),
    })
df = pd.DataFrame(rows)
df.to_csv(os.path.join(TB, 'S4_sensitivity_cohorts.csv'), index=False, encoding='utf-8-sig')

# markdown
fmt = df.copy()
for c in ['Logistic AUC', 'XGBoost AUC', 'Random forest AUC', 'Clinical AUC']:
    fmt[c] = fmt[c].map(lambda v: '%.4f' % v)
lines = []
lines.append('### Table S4. Sensitivity analysis: cohort composition and T1c acquisition geometry')
lines.append('')
lines.append('Same bias-resistant pipeline (identical seed and CV structure) as the primary analysis. '
             'Date-anomaly cohorts of the previous revision are no longer needed: the four affected cases '
             'were resolved (2 re-processed with correct pre-operative series, 1 surgery date corrected, '
             '1 removed from the clinical database).')
lines.append('')
lines.append(fmt.to_markdown(index=False))
lines.append('')
lines.append('- Cohort B **is** the primary analysis cohort (n=89); cohort A adds back the single United '
             'Imaging 1.5 T case (P021, MGMT-negative) and documents the effect of its exclusion.')
lines.append('- Cohorts C/D stratify the primary cohort by T1c slice thickness (3D <= 1.2 mm vs 2D >= 5 mm).')
lines.append('- Clinical baseline in this table is estimated from a single repeat of pooled out-of-fold '
             'predictions; hence its absolute value differs marginally from Table 2 (which averages across repeats).')
open(os.path.join(TB, 'S4_sensitivity_cohorts.md'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')

print(df.to_string(index=False))
a = {r['tag']: r['models'] for r in res}
print('\nA vs B deltas:')
for k in ['logistic', 'xgboost', 'randomforest', 'clinical']:
    print('  %-13s %+.4f' % (k, a['B'][k]['auc'] - a['A'][k]['auc']))
print('\nwritten: sensitivity_cohorts.json, S4_sensitivity_cohorts.csv/.md')
