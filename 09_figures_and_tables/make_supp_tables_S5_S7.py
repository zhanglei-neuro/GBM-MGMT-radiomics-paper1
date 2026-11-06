# -*- coding: utf-8 -*-
r"""Paper 1 -- build Supplementary Tables S5 / S6 / S7 (CSV + Markdown)"""
import os, json, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')

ANA = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/analysis'
TAB = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/tables'
os.makedirs(TAB, exist_ok=True)

def save(df, name, title, notes=None):
    df.to_csv(os.path.join(TAB, name + '.csv'), index=False, encoding='utf-8-sig')
    with open(os.path.join(TAB, name + '.md'), 'w', encoding='utf-8') as f:
        f.write('### %s\n\n' % title)
        f.write(df.to_markdown(index=False))
        if notes:
            f.write('\n\n' + notes + '\n')
    print('saved', name)

# ---------------- S5 all baseline variables ----------------
s5 = pd.read_csv(os.path.join(ANA, 'tableS5_baseline_full.csv'))
save(s5, 'S5_baseline_all_variables',
     'Table S5  Between-group tests for all baseline variables (Siemens vs GE)',
     'Continuous variables: median [IQR], Mann-Whitney U test; categorical variables: n (%), chi-square or Fisher exact test. '
     'Every clinical variable in the table (age, sex, KPS, extent of resection, tumour location, tumour side) was complete '
     'for all 89 cases, with no missing values; all 89 cases are IDH wild-type GBM (WHO 2021 definition), so IDH is not a '
     'between-group comparison. TERT status and follow-up measures are outside the scope of this study and were neither '
     'collected nor analysed.'),

# ---------------- S6 GE subgroup diagnostics ----------------
gd = json.load(open(os.path.join(ANA, 'ge_diagnostics.json')))
folds = pd.DataFrame(gd['folds'])
folds = folds.rename(columns={'fold': 'Fold', 'n_train': 'n train', 'n_test': 'n test',
                              'train_pos': 'train MGMT+', 'train_neg': 'train MGMT−',
                              'test_pos': 'test MGMT+', 'test_neg': 'test MGMT−',
                              'n_selected': 'features selected', 'test_auc': 'test AUC'})
folds['test AUC'] = folds['test AUC'].map(lambda x: '%.3f' % x if np.isfinite(x) else 'NA')
save(folds, 'S6a_GE_fold_diagnostics',
     'Table S6a  GE subgroup (n=32, MGMT+ 13 / MGMT- 19), fold-by-fold diagnostics',
     'Pipeline identical to the main analysis (inner univariate |AUC-0.5| ranking -> greedy de-correlation -> at most 40 features -> inner 3-fold tuning). '
     'Each test fold holds only 6-7 cases, of which just 2-3 are MGMT+, so the AUC takes very coarse values.'),

rows = []
if os.path.exists(os.path.join(ANA, 'ge_perm_null.json')):
    gp = json.load(open(os.path.join(ANA, 'ge_perm_null.json')))
    rows.append(['Label-permutation null (100 permutations, same pipeline)', 'observed AUC', '%.4f' % gp['observed_auc']])
    rows.append(['Label-permutation null (100 permutations, same pipeline)', 'null mean +/- SD', '%.4f +/- %.4f' % (gp['null_mean'], gp['null_sd'])])
    rows.append(['Label-permutation null (100 permutations, same pipeline)', 'null 2.5% - 97.5%',
                 '%.4f – %.4f' % (gp['null_p2.5'], gp['null_p97.5'])])
    rows.append(['Label-permutation null (100 permutations, same pipeline)', 'null min - max', '%.4f - %.4f' % (gp['null_min'], gp['null_max'])])
    rows.append(['Label-permutation null (100 permutations, same pipeline)', 'empirical two-sided p', '%.3f' % gp['perm_p_two_sided']])
    rows.append(['Label-permutation null (100 permutations, same pipeline)', 'fraction of null <= observed', '%.3f' % gp['frac_null_below_observed']])
gu = json.load(open(os.path.join(ANA, 'ge_uni_perm.json')))
rows.append(['Permutation test on the median univariable AUC', 'observed median (400 features)', '%.4f' % gu['observed_median_univariate_auc']])
rows.append(['Permutation test on the median univariable AUC', 'null 2.5% - 97.5%',
             '%.4f – %.4f' % (gu['null_median_p2.5'], gu['null_median_p97.5'])])
rows.append(['Permutation test on the median univariable AUC', 'empirical p (1000 permutations)', '%.3f' % gu['perm_p']])
t10 = gd['top10_direction_on_train']
rows.append(['Stability of the inner-screening direction', 'median univariable AUC of the top-10 training-fold features', '%.2f' % t10['median_auc']])
rows.append(['Stability of the inner-screening direction', 'fraction pointing the other way (AUC<0.5)', '%.0f%%' % (100 * t10['pct_below_0.4'])])
s6b = pd.DataFrame(rows, columns=['Analysis', 'Metric', 'Value'])
save(s6b, 'S6b_GE_permutation_null',
     'Table S6b  GE subgroup: label-permutation null and direction stability',
     'Conclusion: the observed GE AUC (0.198) sits in the lower tail of the 100-permutation same-pipeline null '
     '(2.5th percentile 0.253; empirical two-sided p=0.020), so it does **not** lie squarely in the centre of the null. '
     'That null is nevertheless extremely wide (SD 0.134, range 0.162-0.781, and 1 of the 100 permutations already '
     'produced an AUC of 0.198 or below), and each fold has only 2-3 MGMT+ cases in its test set, which means the AUC '
     'this pipeline estimates at n=32 carries almost no information. The univariable direction test (median 0.478, null '
     '0.449-0.559, p=0.386) further shows that there is no stable systematic inverse association. The paper therefore '
     'keeps only the descriptive statement that performance is below chance and the direction unstable: it neither reads '
     'this as a genuine inverse biological association nor claims that the low value is fully explained by chance.'),

# ---------------- S7 LOVO, three models ----------------
lv = json.load(open(os.path.join(ANA, 'lovo_all_models.json')))
rows = []
for tgt, src in [('Siemens', 'GE'), ('GE', 'Siemens')]:
    d = lv[tgt]
    rows.append(['Train %s (n=%d, MGMT+ %d) → Test %s (n=%d, MGMT+ %d)' %
                 (src, d['n_train'], d['train_pos'], tgt, d['n_test'], d['test_pos']),
                 'Logistic regression', '%.3f' % d['logistic_auc'], d['logistic_n_feat']])
    rows.append(['Train %s (n=%d, MGMT+ %d) → Test %s (n=%d, MGMT+ %d)' %
                 (src, d['n_train'], d['train_pos'], tgt, d['n_test'], d['test_pos']),
                 'XGBoost', '%.3f' % d['xgboost_auc'], d['xgboost_n_feat']])
    rows.append(['Train %s (n=%d, MGMT+ %d) → Test %s (n=%d, MGMT+ %d)' %
                 (src, d['n_train'], d['train_pos'], tgt, d['n_test'], d['test_pos']),
                 'Random forest', '%.3f' % d['randomforest_auc'], d['randomforest_n_feat']])
s7 = pd.DataFrame(rows, columns=['Direction', 'Model', 'AUC', 'No. of features selected'])
save(s7, 'S7_LOVO_all_models',
     'Table S7  Leave-one-vendor-out (LOVO) cross-validation, complete results for the three models',
     'Pipeline as in the main analysis (inner screening + inner 3-fold tuning), fitted on all training-vendor data and '
     'tested on the target vendor. In the "train on GE" direction a 40-feature model is fitted to 32 cases (EPV<1), so its '
     'failure cannot be separated from model capacity, vendor effect and sample size; LOVO is therefore reported only as '
     'an exploratory check of cross-scanner consistency, not as external validation.'),

print('\nDONE supplementary tables S5–S7')
