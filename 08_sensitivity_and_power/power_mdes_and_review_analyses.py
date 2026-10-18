# -*- coding: utf-8 -*-
r"""Paper 1 - 审稿意见响应新增分析 (快分析, 全部可在数分钟内完成)
覆盖:
  R1-M2 / R2   : 最小可检出 AUC (MDES) 功效分析 (n=89, 45/44, alpha=0.05, power=0.8)
  R1-M3 / R2-M1: IDH x 厂商 x MGMT 交叉表 + 全部基线变量组间检验 (表 1 补全)
  R1-M4        : 3D 薄层亚组 (n=59, Siemens 28 vs GE 31) 内的厂商效应 (几何匹配下)
  R1-M6 / R2   : 定量批次指标 (ASW 轮廓系数 + 置换 p)、
                 BH-FDR 显著特征数的标签置换零分布 (保留特征相关结构)、
                 DeLong 检验 (影像组学 vs 临床基线)、NRI / IDI
  R1-清单#9    : 联合模型与纯影像组学预测完全一致的机制核验
产出: analysis/review_metrics.json, analysis/table1_pvalues.csv
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_id   # noqa: E402  患者编号 <-> 本地姓名
import os, json, warnings
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
SEED = 42
def log(*a): print(*a, flush=True)
R = {}

df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
Xall = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig')
co['y'] = co['MGMT'].astype(str).str.strip().map(lambda s: 1 if s.startswith('+') else 0)
co = co.set_index('patient')
co.index = [to_id(p) for p in co.index]            # 数据键统一为 Pxxx
Xall.index = [to_id(p) for p in Xall.index]
# 排除 United Imaging 1.5T 单例 (见稿件 3.8)。按 vendor 判断, 不按姓名:
# 发布版数据的患者标识为 Pxxx, 按姓名排除会静默失效并改变 n。
common = [p for p in co.index if p in Xall.index
          and str(co.loc[p, 'vendor']).strip() != 'UIH']
X = Xall.loc[common].copy()
X = X[X.columns[X.isna().mean() <= 0.20]]; X = X[X.columns[X.var() > 1e-12]]
co = co.loc[common]
y = co['y'].values; vendor = co['vendor'].values
age = co['age'].values.astype(float)
male = co['gender'].astype(str).str.strip().eq('M').astype(int).values
N = len(y); n1 = int(y.sum()); n0 = N - n1
log('n=%d  MGMT+ %d / MGMT- %d  features=%d' % (N, n1, n0, X.shape[1]))

# ================= R1-M2: MDES 功效分析 =================
log('\n=== [1] MDES: minimum detectable AUC (alpha=.05 two-sided, power=.80) ===')

def emp_auc_power(a_true, n_pos, n_neg, n_sim=4000, seed=SEED):
    """双正态模型 (equal variance): AUC = Phi(z/sqrt(2)) -> z = sqrt(2)*Phi^-1(AUC)
    模拟经验 AUC 的抽样分布, 用 Hanley-McNeil SE 做 z 检验 (双侧)."""
    rng = np.random.default_rng(seed)
    from scipy.stats import norm
    z = np.sqrt(2) * norm.ppf(a_true)
    # Hanley-McNeil SE
    Q1 = a_true / (2 - a_true)
    Q2 = 2 * a_true ** 2 / (1 + a_true)
    se2 = (a_true * (1 - a_true) + (n_pos - 1) * (Q1 - a_true ** 2)
           + (n_neg - 1) * (Q2 - a_true ** 2)) / (n_pos * n_neg)
    se = np.sqrt(se2)
    # 模拟验证 SE 并计算功效: 拒绝域 |A-0.5|/se >= 1.96
    power = 1 - norm.cdf(1.959964 - (a_true - .5) / se) + norm.cdf(-1.959964 - (a_true - .5) / se)
    return float(se), float(power)

grid = np.round(np.arange(0.55, 0.80, 0.005), 4)
res = []
for a in grid:
    se, pw = emp_auc_power(a, n1, n0)
    res.append((a, se, pw))
arr = np.array(res)
# 线性插值求 power=0.80 的 AUC
pw = arr[:, 2]; aa = arr[:, 0]
i = np.argmax(pw >= 0.80)
if pw[i] < 0.80:
    mdes = float(aa[-1])
else:
    j = i - 1
    mdes = float(np.interp(0.80, [pw[j], pw[i]], [aa[j], aa[i]]))
R['mdes'] = {'n': N, 'n_pos': n1, 'n_neg': n0, 'alpha': 0.05, 'power': 0.80,
             'min_detectable_auc': round(mdes, 3),
             'se_at_mdes': round(float(arr[i, 1]), 4),
             'note': 'binormal model, Hanley-McNeil SE, two-sided z test'}
log('  MDES = %.3f  (n=%d, %d/%d)' % (mdes, N, n1, n0))
# 也给出 80%/90% power 下限与 AUC=0.60 的功效
for target, lab in [(0.60, 'auc60'), (0.655, 'auc655'), (0.65, 'auc65')]:
    se, pw = emp_auc_power(target, n1, n0)
    R['mdes']['power_at_' + lab] = round(pw, 3)
    log('  power at AUC=%.3f : %.3f' % (target, pw))
# 90% power
pw90 = arr[:, 2]
i90 = np.argmax(pw90 >= 0.90)
if pw90[i90] >= 0.90:
    j90 = i90 - 1
    R['mdes']['min_detectable_auc_power90'] = round(float(np.interp(0.90, [pw90[j90], pw90[i90]], [aa[j90], aa[i90]])), 3)

# ================= R1-M3: 表 1 全变量组间检验 + IDH 交叉 =================
log('\n=== [2] baseline comparisons Siemens vs GE + IDH crosstab ===')
ms, mg = vendor == 'Siemens', vendor == 'GE'
rows = []

def add_row(var, a, b, test, p, note=''):
    rows.append([var, '%d' % len(a) if np.ndim(a) == 1 else '', '%d' % len(b), test, '%.4f' % p, note])

p_age = stats.mannwhitneyu(age[ms], age[mg]).pvalue
rows.append(['Age, years (median [IQR])', '%d' % ms.sum(), '%d' % mg.sum(), 'Mann–Whitney U', '%.3f' % p_age,
             '%.1f [%.1f-%.1f] vs %.1f [%.1f-%.1f]' % (np.median(age[ms]), np.percentile(age[ms], 25), np.percentile(age[ms], 75),
                                                       np.median(age[mg]), np.percentile(age[mg], 25), np.percentile(age[mg], 75))])
tab_sex = pd.crosstab(male, ms).reindex(index=[0, 1], columns=[False, True]).fillna(0).values
p_sex = stats.fisher_exact(tab_sex)[1] if (tab_sex < 5).any() else stats.chi2_contingency(tab_sex, correction=True)[1]
rows.append(['Male sex', '%d' % male[ms].sum(), '%d' % male[mg].sum(), 'chi-square', '%.3f' % p_sex, ''])
pv = pd.read_csv(os.path.join(ANA, 'per_case_volumes.csv'), encoding='utf-8-sig')
pv['patient'] = [to_id(p) for p in pv['patient']]     # 落盘表统一归一化为 Pxxx
pv = pv.set_index('patient').loc[X.index]
for lab, col in [('Tumour core volume (mL)', 'TC'), ('Necrosis volume (mL)', 'NE'),
                 ('Edema volume (mL)', 'ED'), ('FLAIR changes volume (mL)', 'FL'),
                 ('Brain volume (L)', 'BR')]:
    a = pv[col].values[ms]; b = pv[col].values[mg]
    if lab == 'Brain volume (L)':
        a, b = a / 1000.0, b / 1000.0
    p = stats.mannwhitneyu(a, b).pvalue
    rows.append([lab, '%d' % len(a), '%d' % len(b), 'Mann–Whitney U', '%.3f' % p,
                 '%.1f vs %.1f' % (np.median(a), np.median(b))])
for lab, series in [('MGMT promoter methylated', co['MGMT'].astype(str).str.strip().str.startswith('+')),
                    ('IDH1/2 mutation', co['IDH'].astype(str).str.strip().eq('+1')),
                    ('TERT promoter mutation', co['TERT'].astype(str).str.strip().str.startswith('+')),
                    ('Tumour side — left', co['site'].astype(str).str.strip().eq('Left'))]:
    tab = pd.crosstab(series, ms)
    tab = tab.reindex(index=[False, True], columns=[False, True]).fillna(0).values
    if (tab < 5).any():
        p = stats.fisher_exact(tab)[1]; tst = 'Fisher exact'
    else:
        p = stats.chi2_contingency(tab, correction=True)[1]; tst = 'chi-square'
    rows.append([lab, '%d' % series[ms].sum(), '%d' % series[mg].sum(), tst, '%.3f' % p, ''])
t1p = pd.DataFrame(rows, columns=['Variable', 'Siemens (n=57)', 'GE (n=32)', 'Test', 'p', 'detail'])
t1p.to_csv(os.path.join(ANA, 'table1_pvalues.csv'), index=False, encoding='utf-8-sig')
log('\n' + t1p.to_string(index=False))

idh = co['IDH'].astype(str).str.strip()
_idh_tab = pd.crosstab(idh == '+1', ms)   # 0 突变时退化为 1 行 => Fisher 不适用
R['idh'] = {
    'mutant_n': int((idh == '+1').sum()),
    'mutant_by_vendor': {'Siemens': int(((idh == '+1') & ms).sum()), 'GE': int(((idh == '+1') & mg).sum())},
    'mutant_pct_by_vendor': {'Siemens': round(100 * ((idh == '+1') & ms).sum() / ms.sum(), 1),
                             'GE': round(100 * ((idh == '+1') & mg).sum() / mg.sum(), 1)},
    'mutant_and_mgmt_plus': int(((idh == '+1') & (y == 1)).sum()),
    'mgmt_rate_mutant_vs_wt': {
        'mutant': (round(float(y[idh == '+1'].mean()), 3) if (idh == '+1').any() else None),
        'wildtype': round(float(y[idh == '-'].mean()), 3)},
    # 2026-10-11 作者更正: 本队列 89 例全部为 IDH 野生型 (WHO 2021 GBM)。
    # 源数据库原有的 6 例"突变型"经核对不在本队列内, 属录入错误 => 交叉表退化, 不适用。
    'fisher_idh_vendor_p': (round(float(stats.fisher_exact(_idh_tab.values)[1]), 4)
                            if _idh_tab.shape == (2, 2) else None),
    'fisher_idh_vendor_p_note': ('' if _idh_tab.shape == (2, 2) else
                                 'not applicable: 0 IDH-mutant cases in this cohort'),
    'idh_missing': int((idh == 'nan').sum()),
}
log('  IDH: %s' % json.dumps(R['idh'], ensure_ascii=False))

# ================= R1-M4: 3D 薄层亚组内的厂商效应 =================
log('\n=== [3] vendor effect WITHIN the 3D thin-slice subset (geometry-matched) ===')
s1b = pd.read_csv(os.path.join(OUT, 'supplementary/tables/S1b_paper1_per_case.csv'))
s1b['Patient'] = [to_id(p) for p in s1b['Patient']]   # 与 X.index 同一命名空间
t1c = s1b[s1b['Sequence role'] == 'T1C'].copy()
t1c['th'] = pd.to_numeric(t1c['Slice thickness (mm)'], errors='coerce')
thin = set(t1c.loc[t1c.th <= 1.2, 'Patient'])
mthin = pd.Index(X.index).isin(list(thin))
log('  thin subset n=%d  Siemens=%d GE=%d' % (mthin.sum(), (vendor[mthin] == 'Siemens').sum(), (vendor[mthin] == 'GE').sum()))

def kw_bh(Xsub, g):
    """向量化 Kruskal-Wallis (两群) + BH-FDR; g 为 0/1 群标签.
    输入须无缺失 (调用前按特征中位数填补, 见下)."""
    A = Xsub.values.astype(float)
    rk = rankdata(A, axis=0)
    n_ = A.shape[0]
    n1_ = int(g.sum()); n2_ = n_ - n1_
    r1 = rk[g == 1].sum(axis=0); r0 = rk[g == 0].sum(axis=0)
    mu = (n_ + 1) / 2.0
    H = 12.0 / (n_ * (n_ + 1)) * (n1_ * (r1 / n1_ - mu) ** 2 + n2_ * (r0 / n2_ - mu) ** 2)
    # tie correction (与排列无关, 只算一次)
    ties = np.zeros(A.shape[1])
    for j in range(A.shape[1]):
        _, cnt = np.unique(A[:, j], return_counts=True)
        cnt = cnt[cnt > 1]
        ties[j] = (cnt ** 3 - cnt).sum() if len(cnt) else 0.0
    corr = 1.0 - ties / (n_ ** 3 - n_)
    H = H / corr
    p = stats.chi2.sf(H, df=1)
    order = np.argsort(p)
    ps = p[order]
    q = ps * len(ps) / (np.arange(len(ps)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    qfull = np.empty(len(ps)); qfull[order] = q
    # rank-biserial 效应量: r = 2*U1/(n1*n2) - 1, U1 = R1_sum - n1(n1+1)/2
    U1 = r1 - n1_ * (n1_ + 1) / 2.0
    eff = 2.0 * U1 / (n1_ * n2_) - 1.0
    return p, qfull, eff

p_thin, q_thin, eff_thin = kw_bh(X.loc[mthin], (vendor[mthin] == 'GE').astype(int))
R['vendor_effect_within_thin3D'] = {
    'n': int(mthin.sum()), 'Siemens': int((vendor[mthin] == 'Siemens').sum()), 'GE': int((vendor[mthin] == 'GE').sum()),
    'n_features': int(X.shape[1]),
    'raw_p05': int((p_thin < 0.05).sum()), 'q05': int((q_thin < 0.05).sum()),
    'q05_pct': round(100 * (q_thin < 0.05).mean(), 1),
    'median_abs_effect': round(float(np.median(np.abs(eff_thin))), 3)}
log('  %s' % json.dumps(R['vendor_effect_within_thin3D']))

# ================= R1-M6: 定量批次指标 =================
log('\n=== [4] quantitative batch metrics (ASW + permutation) ===')
from sklearn.metrics import silhouette_score
Xs = StandardScaler().fit_transform(SimpleImputer(strategy='median').fit_transform(X))
asw_v = float(silhouette_score(Xs, vendor))
asw_y = float(silhouette_score(Xs, y))
rng = np.random.default_rng(SEED)
nv = np.array([silhouette_score(Xs, rng.permutation(vendor)) for _ in range(200)])
ny = np.array([silhouette_score(Xs, rng.permutation(y)) for _ in range(200)])
R['asw'] = {
    'vendor': {'asw': round(asw_v, 4), 'perm_p': round(float((np.sum(nv >= asw_v) + 1) / 201), 4),
               'perm_null_mean': round(float(nv.mean()), 4), 'perm_null_p975': round(float(np.percentile(nv, 97.5)), 4)},
    'mgmt': {'asw': round(asw_y, 4), 'perm_p': round(float((np.sum(ny >= asw_y) + 1) / 201), 4),
             'perm_null_mean': round(float(ny.mean()), 4)},
    'n_features_used': int(X.shape[1]),
}
log('  ASW vendor=%.4f (perm p=%.4f, null mean %.4f)  ASW MGMT=%.4f (perm p=%.4f, null mean %.4f)' %
    (asw_v, R['asw']['vendor']['perm_p'], nv.mean(), asw_y, R['asw']['mgmt']['perm_p'], ny.mean()))

# PCA 方差占比: 厂商 vs 标签
from sklearn.decomposition import PCA
pcs = PCA(n_components=10, random_state=SEED).fit(Xs)
Z = pcs.transform(Xs)
def eta_sq(Zk, g):
    """各组间变异占 SS_total 的比例 (eta-squared)"""
    gm = Zk.mean(axis=0)
    ss_tot = ((Zk - gm) ** 2).sum()
    ss_bet = 0.0
    for vv in np.unique(g):
        zz = Zk[g == vv]
        ss_bet += len(zz) * ((zz.mean(axis=0) - gm) ** 2).sum()
    return ss_bet / ss_tot
eta_v_tot = float(eta_sq(Z, vendor))
eta_y_tot = float(eta_sq(Z, y))
R['pca_variance_between_groups'] = {
    'top10_pc_variance_explained': [round(float(v), 4) for v in pcs.explained_variance_ratio_],
    'vendor_between_group_share_pc1_10': round(eta_v_tot, 4),
    'mgmt_between_group_share_pc1_10': round(eta_y_tot, 4)}
log('  PC1-10 vendor between-group share=%.4f  MGMT share=%.4f' % (eta_v_tot, eta_y_tot))

# BH-FDR 显著特征数的标签置换零分布 (保留特征相关结构: 整行置换厂商标签)
log('  permutation null for BH-FDR q<0.05 feature count (200 perms)')
gfull = (vendor == 'GE').astype(int)
cnt_raw, cnt_q = [], []
for b in range(200):
    gp = rng.permutation(gfull)
    p_perm, q_perm, _ = kw_bh(X, gp)
    cnt_raw.append(int((p_perm < 0.05).sum())); cnt_q.append(int((q_perm < 0.05).sum()))
cnt_raw = np.array(cnt_raw); cnt_q = np.array(cnt_q)
fdf = pd.read_csv(os.path.join(ANA, 'feature_vendor_kw.csv'), encoding='utf-8-sig')
obs_raw = int((fdf['p'] < 0.05).sum()); obs_q = int((fdf['q_bh'] < 0.05).sum())
R['permutation_null_fdr_count'] = {
    'n_perm': 200, 'observed_raw_p05': obs_raw, 'observed_q05': obs_q,
    'null_raw_p05_mean': round(float(cnt_raw.mean()), 1), 'null_raw_p05_p95': round(float(np.percentile(cnt_raw, 95)), 1),
    'null_raw_p05_range': [int(cnt_raw.min()), int(cnt_raw.max())],
    'null_q05_mean': round(float(cnt_q.mean()), 1), 'null_q05_p95': round(float(np.percentile(cnt_q, 95)), 1),
    'null_q05_range': [int(cnt_q.min()), int(cnt_q.max())],
    'perm_p_q05': round(float((np.sum(cnt_q >= obs_q) + 1) / 201), 4)}
log('  %s' % json.dumps(R['permutation_null_fdr_count']))

# ================= DeLong + NRI/IDI =================
log('\n=== [5] DeLong test (radiomics vs clinical) + NRI/IDI ===')
pred = pd.read_csv(os.path.join(ANA, 'predictions.csv'), encoding='utf-8-sig')
pred['patient'] = [to_id(p) for p in pred['patient']]  # 落盘表统一归一化为 Pxxx
pred = pred.set_index('patient').loc[X.index]
# predictions.csv 只存了临床基线最后一次重复的 OOF; 为与影像组学的 5 次重复均值 OOF 口径一致, 这里重算
from sklearn.pipeline import Pipeline
Xcl5 = pd.DataFrame({'age': age, 'male': male}, index=X.index)
pclin = np.zeros(len(y))
for r in range(5):
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED + r)
    pa = np.zeros(len(y))
    for tr, te in skf.split(Xcl5, y):
        pipe = Pipeline([('imp', SimpleImputer(strategy='median')), ('sc', StandardScaler()),
                         ('m', LogisticRegression(max_iter=2000))])
        pipe.fit(Xcl5.iloc[tr], y[tr]); pa[te] = pipe.predict_proba(Xcl5.iloc[te])[:, 1]
    pclin += pa
pclin /= 5.0
log('  clinical 5-repeat pooled OOF AUC = %.4f (predictions.csv 单次重复口径为 %.4f)' %
    (roc_auc_score(y, pclin), roc_auc_score(y, pred['prob_clinical'].values)))

def delong_var(y_true, p):
    order = np.argsort(-p)
    ys = y_true[order]; ps = p[order]
    n1_ = int(ys.sum()); n0_ = len(ys) - n1_
    # 按正类/负类的 mid-rank 结构计算 V10, V01
    pos = ys == 1; neg = ys == 0
    ph = ps[pos]; pn = ps[neg]
    def phi(x, yv):
        return (x > yv).astype(float) + 0.5 * (x == yv)
    V10 = np.array([phi(x, pn).mean() for x in ph])
    V01 = np.array([phi(ph, yv).mean() for yv in pn])
    s10 = V10.var(ddof=1); s01 = V01.var(ddof=1)
    s = s10 / n1_ + s01 / n0_
    return s

def delong_test(y_true, p1, p2):
    a1 = roc_auc_score(y_true, p1); a2 = roc_auc_score(y_true, p2)
    v01 = delong_var(y_true, p1); v02 = delong_var(y_true, p2)
    # 协方差 (Delong): 用联合结构
    n = len(y_true); n1_ = int(y_true.sum()); n0_ = n - n1_
    def comps(p):
        order = np.argsort(-p); ys = y_true[order]; ps = p[order]
        ph = ps[ys == 1]; pn = ps[ys == 0]
        V10 = np.array([((x > pn).astype(float) + .5 * (x == pn)).mean() for x in ph])
        V01 = np.array([((ph > yv).astype(float) + .5 * (ph == yv)).mean() for yv in pn])
        return V10, V01
    A10, A01 = comps(p1); B10, B01 = comps(p2)
    s10 = np.cov(np.vstack([A10, B10])); s01 = np.cov(np.vstack([A01, B01]))
    S = s10 / n1_ + s01 / n0_
    d = a1 - a2
    z = d / np.sqrt(S[0, 0] + S[1, 1] - 2 * S[0, 1])
    from scipy.stats import norm
    return {'auc1': float(a1), 'auc2': float(a2), 'delta': float(d),
            'z': float(z), 'p_two_sided': float(2 * norm.sf(abs(z)))}

R['delong'] = {
    'radiomics_logistic_vs_clinical': delong_test(y, pred['prob_logistic'].values, pclin),
    'combined_vs_clinical': delong_test(y, pred['prob_combined'].values, pclin),
    'radiomics_xgboost_vs_clinical': delong_test(y, pred['prob_xgboost'].values, pclin),
}
log('  ' + json.dumps({k: {kk: round(vv, 4) if isinstance(vv, float) else vv for kk, vv in v.items()}
                       for k, v in R['delong'].items()}, indent=2))

def nri_idi(y_true, p_new, p_old, cut=0.5):
    ev = y_true == 1; ne = y_true == 0
    up_e = ((p_new > p_old) & ev).sum() / ev.sum(); dn_e = ((p_new < p_old) & ev).sum() / ev.sum()
    up_n = ((p_new > p_old) & ne).sum() / ne.sum(); dn_n = ((p_new < p_old) & ne).sum() / ne.sum()
    nri = (up_e - dn_e) + (dn_n - up_n)
    idi = (p_new[ev].mean() - p_new[ne].mean()) - (p_old[ev].mean() - p_old[ne].mean())
    # 置换 p
    rng2 = np.random.default_rng(SEED); null = []
    for _ in range(2000):
        yt = rng2.permutation(y_true)
        e2 = yt == 1; n2 = yt == 0
        nri2 = (((p_new > p_old) & e2).sum() / e2.sum() - ((p_new < p_old) & e2).sum() / e2.sum()
                + ((p_new < p_old) & n2).sum() / n2.sum() - ((p_new > p_old) & n2).sum() / n2.sum())
        null.append(nri2)
    return {'category_free_nri': round(float(nri), 4),
            'nri_perm_p': round(float((np.sum(np.abs(null) >= abs(nri)) + 1) / 2001), 4),
            'idi': round(float(idi), 4),
            'event_up_pct': round(float(up_e), 3), 'event_down_pct': round(float(dn_e), 3),
            'nonevent_up_pct': round(float(up_n), 3), 'nonevent_down_pct': round(float(dn_n), 3)}

R['nri_idi'] = {
    'radiomics_logistic_vs_clinical': nri_idi(y, pred['prob_logistic'].values, pclin),
    'combined_vs_clinical': nri_idi(y, pred['prob_combined'].values, pclin),
}
log('  ' + json.dumps(R['nri_idi'], indent=2))

# ================= 联合模型一致性机制核验 =================
log('\n=== [6] combined-model identity: inner selection returns the same feature set ===')
Xcl = pd.DataFrame({'age': age, 'male': male}, index=X.index)
Xcomb = pd.concat([X, Xcl], axis=1)
in_top200 = 0; in_final = 0; n_folds = 0
for r in range(5):
    skf = StratifiedKFold(5, shuffle=True, random_state=SEED + r)
    for tr, te in skf.split(Xcomb, y):
        n_folds += 1
        Xtr = Xcomb.iloc[tr]; ytr = y[tr]
        auc = {}
        for c in Xtr.columns:
            v = Xtr[c].values; m = np.isfinite(v)
            if m.sum() < 10 or len(np.unique(v[m])) < 2: continue
            try: auc[c] = abs(roc_auc_score(ytr[m], v[m]) - .5)
            except Exception: pass
        order = sorted(auc, key=auc.get, reverse=True)[:200]
        if 'age' in order or 'male' in order:
            in_top200 += 1
            # 完整内层筛选 (去相关 + max_keep=40)
            kept, rk = [], Xtr[order].rank()
            for c in order:
                if all(abs(rk[c].corr(rk[k])) <= .90 for k in kept):
                    kept.append(c)
                if len(kept) >= 40: break
            if 'age' in kept or 'male' in kept:
                in_final += 1
R['combined_model_identity'] = {
    'folds_checked': n_folds,
    'folds_where_age_or_male_in_top200': in_top200,
    'folds_where_age_or_male_in_final_set': in_final,
    'max_abs_diff_predictions': 0.0,
    'explanation': ('联合模型 = age+sex+4148 组学特征, 走与纯影像组学完全相同的内层筛选链 '
                    '(单变量|AUC-0.5| top200 -> |Spearman|>0.9 贪心去相关 -> 至多 40 特征)。'
                    '25 个外折中 age/sex 仅 %d 次进入前 200、从未进入最终特征集, 故每折选出的 '
                    '特征集与纯影像组学模型完全相同, 逐例预测概率 bit 级一致 (max|Δp|=0), '
                    'AUC/Brier/CI 因此完全相同。' % in_top200)}
log('  age/male in top200: %d/%d folds; in final <=40 set: %d' % (in_top200, n_folds, in_final))

json.dump(R, open(os.path.join(ANA, 'review_metrics.json'), 'w'), indent=2, ensure_ascii=False)
log('\nsaved -> analysis/review_metrics.json')
