# -*- coding: utf-8 -*-
r"""Paper 1 -- 生成补充表 S5 / S6 / S7 (CSV + Markdown)"""
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

# ---------------- S5 基线全变量 ----------------
s5 = pd.read_csv(os.path.join(ANA, 'tableS5_baseline_full.csv'))
save(s5, 'S5_baseline_all_variables',
     '表 S5  队列基线特征全变量组间检验（Siemens vs GE）',
     '连续变量：中位数 [四分位距]，Mann–Whitney U 检验；分类变量：n (%)，χ² 或 Fisher 精确检验。'
     '表内全部临床变量（年龄、性别、KPS、切除程度、肿瘤部位、肿瘤侧别）在 89 例中均完整记录、'
     '不存在缺失；全队列 89 例均为 IDH 野生型 GBM（WHO 2021 定义实体），故 IDH 不构成组间比较项。'
     'TERT 状态与随访指标不属于本研究范围，未纳入采集与分析。')

# ---------------- S6 GE 亚组诊断 ----------------
gd = json.load(open(os.path.join(ANA, 'ge_diagnostics.json')))
folds = pd.DataFrame(gd['folds'])
folds = folds.rename(columns={'fold': 'Fold', 'n_train': 'n train', 'n_test': 'n test',
                              'train_pos': 'train MGMT+', 'train_neg': 'train MGMT−',
                              'test_pos': 'test MGMT+', 'test_neg': 'test MGMT−',
                              'n_selected': 'features selected', 'test_auc': 'test AUC'})
folds['test AUC'] = folds['test AUC'].map(lambda x: '%.3f' % x if np.isfinite(x) else 'NA')
save(folds, 'S6a_GE_fold_diagnostics',
     '表 S6a  GE 亚组（n=32，MGMT+ 13 / MGMT− 19）逐折诊断',
     '管线与主分析完全一致（内层单变量 |AUC−0.5| 排序 → 贪心去相关 → 至多 40 特征 → 内层 3 折调参）。'
     '测试折仅含 6–7 例、其中 MGMT+ 仅 2–3 例，AUC 取值格点极粗。')

rows = []
if os.path.exists(os.path.join(ANA, 'ge_perm_null.json')):
    gp = json.load(open(os.path.join(ANA, 'ge_perm_null.json')))
    rows.append(['标签置换零分布（100 次，同管线）', 'observed AUC', '%.4f' % gp['observed_auc']])
    rows.append(['标签置换零分布（100 次，同管线）', 'null mean ± SD', '%.4f ± %.4f' % (gp['null_mean'], gp['null_sd'])])
    rows.append(['标签置换零分布（100 次，同管线）', 'null 2.5% – 97.5%',
                 '%.4f – %.4f' % (gp['null_p2.5'], gp['null_p97.5'])])
    rows.append(['标签置换零分布（100 次，同管线）', 'null min – max', '%.4f – %.4f' % (gp['null_min'], gp['null_max'])])
    rows.append(['标签置换零分布（100 次，同管线）', '经验双侧 p', '%.3f' % gp['perm_p_two_sided']])
    rows.append(['标签置换零分布（100 次，同管线）', 'null ≤ observed 的比例', '%.3f' % gp['frac_null_below_observed']])
gu = json.load(open(os.path.join(ANA, 'ge_uni_perm.json')))
rows.append(['单变量 AUC 中位数置换检验', '观察中位数（400 特征）', '%.4f' % gu['observed_median_univariate_auc']])
rows.append(['单变量 AUC 中位数置换检验', '零分布 2.5% – 97.5%',
             '%.4f – %.4f' % (gu['null_median_p2.5'], gu['null_median_p97.5'])])
rows.append(['单变量 AUC 中位数置换检验', '经验 p（1000 次）', '%.3f' % gu['perm_p']])
t10 = gd['top10_direction_on_train']
rows.append(['内层筛选方向稳定性', '训练折 top10 特征中位单变量 AUC', '%.2f' % t10['median_auc']])
rows.append(['内层筛选方向稳定性', '其中方向相反（AUC<0.5）的比例', '%.0f%%' % (100 * t10['pct_below_0.4'])])
s6b = pd.DataFrame(rows, columns=['分析', '指标', '数值'])
save(s6b, 'S6b_GE_permutation_null',
     '表 S6b  GE 亚组标签置换零分布与方向稳定性',
     '结论：GE 亚组的观察 AUC（0.198）位于 100 次同管线标签置换零分布的下尾（2.5% 分位 0.253；'
     '经验双侧 p=0.020），**并非完全落在零分布中心区**；但该零分布极宽（SD 0.134，范围 0.162–0.781，'
     '100 次置换中已有 1 次产生 ≤0.198 的 AUC），且各折测试集仅 2–3 例 MGMT+，说明该管线的 AUC 估计'
     '在 n=32 下几乎不携带信息。单变量层面的方向置换检验（中位数 0.478，零分布 0.449–0.559，p=0.386）'
     '进一步表明不存在稳定的系统性反向关联。故本文仅保留"性能低于随机且方向不稳定"这一描述性表述，'
     '既不解释为真实的反向生物学关联，也不声称该低值完全由随机波动解释。')

# ---------------- S7 LOVO 三模型 ----------------
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
     '表 S7  留一厂商交叉验证（LOVO）三模型完整结果',
     '管线与主分析一致（内层筛选 + 内层 3 折调参），在训练厂商全量数据上拟合后在目标厂商上测试。'
     '「训练 GE」方向下为 32 例拟合 40 维模型（EPV<1），其失败与模型容量、厂商效应、样本量不可分离，'
     '故 LOVO 仅作为探索性跨设备一致性参考，不作为外部验证。')

print('\nDONE supplementary tables S5–S7')
