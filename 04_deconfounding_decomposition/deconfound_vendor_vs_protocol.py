# -*- coding: utf-8 -*-
r"""Paper 1 -- R1-M4(i): "厂商效应"在调整采集参数后的剩余效应

主分析（与稿件 §3.x / 图 S4 / 答复函 R1-M4 完全一致的注册口径）
--------------------------------------------------------------
对每个组学特征 f, 以其来源序列的采集参数为协变量, 拟合

    z(f) ~ vendor + type2D + slice_thickness + TR + TE + TI

取 vendor 系数的 t 检验 p 值, 与"仅 vendor"的 Kruskal-Wallis p 值对比,
统计 BH-FDR q<0.05 的特征数变化, 从而量化"厂商"在扣除采集参数后的剩余效应。
预期结果（稿件引用值）：3382 (81.5%) -> 259 (6.2%)；配对 224 / 3158 (93.4%)；
偏 R² 中位数 vendor|acquisition 0.016 vs acquisition|vendor 0.098（约 6 倍）。

敏感性分析（不作为稿件主数字，仅存档于 sensitivity_inplane_covariate）
----------------------------------------------------------------------
在上述协变量集基础上追加**面内像素间距**
（= (Pixel spacing row + Pixel spacing col)/2, 取自该特征实际所在序列的逐例
DICOM 实测参数, 补充表 S1b）。纳入理由: 面内间距在本队列为 0.34-1.30 mm,
与厂商相关 (r=0.375), 遗漏该协变量可能影响"厂商"的净效应估计。
该口径下 adjusted q<0.05 = 362 (8.7%)，结论方向不变（约九成厂商效应仍可由
采集参数解释），故稿件正文仍报告 5 协变量注册口径。

方差分解: 对全特征的 median partial R² (vendor | acquisition) 与
(acquisition | vendor)（主口径），以及 (in-plane | 其余)（敏感性口径）。

输出: analysis/m4_vendor_vs_protocol.json
"""
import os, json, warnings
import numpy as np, pandas as pd
from scipy import stats
import statsmodels.api as sm
warnings.filterwarnings('ignore')

SRC = r'C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts'
OUT = r'D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1'
ANA = os.path.join(OUT, 'analysis')
TAB = os.path.join(OUT, 'supplementary', 'tables')
SEED = 42
rng = np.random.default_rng(SEED)

COV_BASE = ['type2D', 'thk', 'TR', 'TE', 'TI']          # 稿件注册口径（5 协变量）
COV_SENS = COV_BASE + ['inplane']                       # 敏感性口径（追加面内像素间距）


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


# ---------- 特征矩阵 ----------
df = pd.read_csv(os.path.join(SRC, 'radiomics_features_long.csv'), encoding='utf-8-sig')
df['feat'] = df['roi'] + '|' + df['seq'] + '|' + df['feature']
X = df.pivot_table(index='patient', columns='feat', values='value', aggfunc='first')
X = X[X.columns[X.isna().mean() <= 0.20]]
X = X[X.columns[X.var() > 1e-12]]
X = X.loc[:, X.notna().any()]

co = pd.read_csv(os.path.join(SRC, 'cohort_MGMT.csv'), encoding='utf-8-sig').set_index('patient')
# 排除 United Imaging 单例 (n=1): 无法作平台层面推断, 且 ComBat 单例 batch 无法估计批参数。
# 按厂商排除而非按姓名, 以免因姓名脱敏 (<PATIENT>) 失效而误纳入 => 保证 n=89 与稿件一致。
common = [p for p in co.index if p in X.index and str(co.loc[p, 'vendor']).strip() != 'UIH']
X = X.loc[common]
y_vendor = (co.loc[common, 'vendor'] == 'Siemens').astype(int)   # 1 = Siemens
print('cohort after excluding UIH: n=%d' % len(common))

# ---------- 逐例采集参数 (取实际用于特征提取的那条序列) ----------
acq = pd.read_csv(os.path.join(TAB, 'S1b_paper1_per_case.csv'), encoding='utf-8-sig')
acq.columns = [c.strip() for c in acq.columns]
acq = acq.rename(columns={'Sequence role': 'seq', 'Acquisition type': 'atype'})
acq['Patient'] = acq['Patient'].astype(str).str.strip()
VCOLS = ['TR (ms)', 'TE (ms)', 'TI (ms)', 'Slice thickness (mm)',
         'Pixel spacing row (mm)', 'Pixel spacing col (mm)']
for c in VCOLS:
    acq[c] = pd.to_numeric(acq[c], errors='coerce')
acq['type2D'] = acq['atype'].astype(str).str.contains('2D').astype(float)
# 面内像素间距 = 行/列间距均值
acq['inplane'] = acq[['Pixel spacing row (mm)', 'Pixel spacing col (mm)']].mean(axis=1)

# 每个 (patient, seq) 一行
A = acq.set_index(['Patient', 'seq'])
print('acquisition rows:', len(acq), 'patients:', acq.Patient.nunique())
print(acq.groupby('seq')[VCOLS + ['type2D']].median().to_string())

# ---------- 设计矩阵 (逐特征换参数) ----------
SEQ_OF = {}
for c in X.columns:
    parts = c.split('|')
    s = parts[1]
    SEQ_OF[c] = 'T1C' if s == 'T1CE' else s      # 特征侧命名 T1CE, 采集表侧命名 T1C


def design_for(seq, idx):
    """返回该序列对应的协变量矩阵 (n x k)"""
    sub = A.xs(seq, level='seq') if seq in set(acq['seq']) else None
    if sub is None or len(sub) == 0:
        return None
    sub = sub.reindex(idx)
    D = pd.DataFrame(index=idx)
    D['type2D'] = sub['type2D'].values
    D['thk'] = sub['Slice thickness (mm)'].values
    D['TR'] = sub['TR (ms)'].values
    D['TE'] = sub['TE (ms)'].values
    D['TI'] = sub['TI (ms)'].values
    D['inplane'] = sub['inplane'].values
    return D


# 检查序列名是否匹配
avail_roles = set(acq['seq'].astype(str))
need_roles = set(SEQ_OF.values())
print('feature seq roles:', sorted(need_roles), '| acquisition roles:', sorted(avail_roles))
missing_roles = need_roles - avail_roles
if missing_roles:
    print('WARNING missing roles:', missing_roles)


# 标准化协变量 (连续列)
def zs(s):
    s = pd.Series(s, dtype=float)
    sd = s.std(ddof=0)
    return (s - s.mean()) / (sd if sd > 0 else 1.0)


idx = X.index
cov_cache = {}
for seq in need_roles:
    D = design_for(seq, idx)
    if D is None:
        continue
    for c in ['thk', 'TR', 'TE', 'TI', 'inplane']:
        D[c] = zs(D[c])
    cov_cache[seq] = D

Xv = X.values.astype(float)


def bh_run(covs):
    """给定协变量集合, 返回 (adj_p, pr2_vendor, pr2_acq, pr2_inplane)。
    pr2_* 仅当相应列在模型中时才有意义; inplane 不在 covs 中时返回全 NaN。"""
    adj_p = np.full(X.shape[1], np.nan)
    pr2_vendor = np.full(X.shape[1], np.nan)   # vendor | acquisition
    pr2_acq = np.full(X.shape[1], np.nan)      # acquisition | vendor
    pr2_inplane = np.full(X.shape[1], np.nan)  # in-plane | vendor + 其余采集参数
    has_ip = 'inplane' in covs

    def arsq(res, res_r):
        if res.rsquared is None or res_r.rsquared is None:
            return np.nan
        return (res.rsquared - res_r.rsquared) / (1 - res_r.rsquared) \
            if 1 - res_r.rsquared > 1e-9 else np.nan

    for j, col in enumerate(X.columns):
        seq = SEQ_OF[col]
        D = cov_cache.get(seq)
        if D is None:
            continue
        v = Xv[:, j]
        m = np.isfinite(v) & np.isfinite(D[covs].values).all(axis=1)
        if m.sum() < 30:
            continue
        try:
            yy = zs(pd.Series(v[m]))
            Dm = D.loc[m]
            base = Dm[covs].copy()
            full = base.copy(); full['vendor'] = y_vendor.loc[m].values
            rf = sm.OLS(yy.values, sm.add_constant(full.astype(float))).fit()
            rr = sm.OLS(yy.values, sm.add_constant(base.astype(float))).fit()
            adj_p[j] = rf.pvalues['vendor']
            pr2_vendor[j] = arsq(rf, rr)
            # acquisition | vendor
            rr2 = sm.OLS(yy.values,
                         sm.add_constant(pd.DataFrame(
                             {'vendor': y_vendor.loc[m].values}).astype(float))).fit()
            pr2_acq[j] = arsq(rf, rr2)
            if has_ip:
                ri = sm.OLS(yy.values,
                            sm.add_constant(full.drop(columns=['inplane']).astype(float))).fit()
                pr2_inplane[j] = arsq(rf, ri)
        except Exception:
            pass
    return adj_p, pr2_vendor, pr2_acq, pr2_inplane


# ---------- 仅 vendor 的 KW p 值 (基线) ----------
kw_p = np.full(X.shape[1], np.nan)
for j in range(X.shape[1]):
    v = Xv[:, j]; m = np.isfinite(v)
    if m.sum() < 10:
        continue
    a = v[m & (y_vendor.values == 1)]; b = v[m & (y_vendor.values == 0)]
    if len(a) < 5 or len(b) < 5:
        continue
    kw_p[j] = stats.kruskal(a, b).pvalue

kw_q = bh(np.nan_to_num(kw_p, nan=1.0))

results = {'n': int(len(idx)), 'n_features': int(X.shape[1]),
           'note': ('primary model (registered, matches manuscript & Figure S4): '
                    'z(feature) ~ vendor + type2D + thickness + TR + TE + TI')}

# ---------- 主分析: 5 协变量注册口径 ----------
adj_p, pr2_vendor, pr2_acq, _ = bh_run(COV_BASE)
adj_q = bh(np.nan_to_num(adj_p, nan=1.0))
valid = np.isfinite(kw_p) & np.isfinite(adj_p)

results['unadjusted_vendor_only'] = {
    'raw_p05': int(np.nansum(kw_p < 0.05)),
    'q05': int(np.nansum(kw_q < 0.05)),
    'pct_q05': float(100 * np.nansum(kw_q < 0.05) / len(kw_q)),
}
results['adjusted_for_acquisition'] = {
    'raw_p05': int(np.nansum(adj_p < 0.05)),
    'q05': int(np.nansum(adj_q < 0.05)),
    'pct_q05': float(100 * np.nansum(adj_q < 0.05) / len(adj_q)),
}
results['paired_on_valid'] = {
    'n_valid_features': int(valid.sum()),
    'q05_both': int(np.sum((kw_q[valid] < 0.05) & (adj_q[valid] < 0.05))),
    'q05_unadj_only': int(np.sum((kw_q[valid] < 0.05) & (adj_q[valid] >= 0.05))),
    'q05_adj_only': int(np.sum((kw_q[valid] >= 0.05) & (adj_q[valid] < 0.05))),
}
results['partial_r2'] = {
    'vendor_given_acquisition_median': float(np.nanmedian(pr2_vendor)),
    'acquisition_given_vendor_median': float(np.nanmedian(pr2_acq)),
    'vendor_given_acquisition_mean': float(np.nanmean(pr2_vendor)),
    'acquisition_given_vendor_mean': float(np.nanmean(pr2_acq)),
}

# ---------- 敏感性: 追加面内像素间距 (6 协变量) ----------
adj_p6, pr2_vendor6, pr2_acq6, pr2_inplane = bh_run(COV_SENS)
adj_q6 = bh(np.nan_to_num(adj_p6, nan=1.0))
valid6 = np.isfinite(kw_p) & np.isfinite(adj_p6)
results['sensitivity_inplane_covariate'] = {
    'note': ('sensitivity only, NOT the manuscript headline number: '
             'z(feature) ~ vendor + type2D + thickness + TR + TE + TI '
             '+ in-plane pixel spacing'),
    'raw_p05': int(np.nansum(adj_p6 < 0.05)),
    'q05': int(np.nansum(adj_q6 < 0.05)),
    'pct_q05': float(100 * np.nansum(adj_q6 < 0.05) / len(adj_q6)),
    'paired_on_valid': {
        'n_valid_features': int(valid6.sum()),
        'q05_both': int(np.sum((kw_q[valid6] < 0.05) & (adj_q6[valid6] < 0.05))),
        'q05_unadj_only': int(np.sum((kw_q[valid6] < 0.05) & (adj_q6[valid6] >= 0.05))),
        'q05_adj_only': int(np.sum((kw_q[valid6] >= 0.05) & (adj_q6[valid6] < 0.05))),
    },
    'partial_r2': {
        'vendor_given_acquisition_median': float(np.nanmedian(pr2_vendor6)),
        'acquisition_given_vendor_median': float(np.nanmedian(pr2_acq6)),
        'inplane_given_rest_median': float(np.nanmedian(pr2_inplane)),
        'inplane_given_rest_mean': float(np.nanmean(pr2_inplane)),
    },
}

# 共线性诊断 (以 T1C 设计矩阵为例; 两两 Pearson 与是否纳入 inplane 无关)
# 注意: SEQ_OF 已把特征侧 T1CE 映射为采集表侧 T1C, 故 cov_cache 的键是 'T1C'。
D = cov_cache['T1C'] if 'T1C' in cov_cache else cov_cache.get('T1CE')
if D is not None:
    full = D.copy(); full['vendor'] = y_vendor.values
    c = full.corr().abs()
    results['collinearity_max_abs_r_with_vendor'] = \
        {k: float(v) for k, v in c['vendor'].drop('vendor').items()}

json.dump(results, open(os.path.join(ANA, 'm4_vendor_vs_protocol.json'), 'w'), indent=2)
print(json.dumps(results, indent=2))
print('saved -> analysis/m4_vendor_vs_protocol.json')
