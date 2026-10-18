# -*- coding: utf-8 -*-
r"""
批量检查生产序列选择（series_selection.csv）中三个序列是否来自同一次检查。
路径形如 D:\NeuroImagepre\<pat>\<YYYYMMDD>\... 或 D:\NeuroImagepre\<pat>\ImageFileName_xxx\<n>
"""
import pandas as pd, re, os

SEL = r"C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\series_selection.csv"
s = pd.read_csv(SEL)
s = s[s.patient.notna()]


def get_date(p):
    p = str(p)
    m = re.search(r"\\(\d{8})\\", p)
    return m.group(1) if m else "<no-date-layer>"


rows = []
for _, r in s.iterrows():
    d = {c: get_date(r[c]) for c in ["t1ce", "flair", "t1"]}
    uniq = set(d.values())
    rows.append({
        "patient": r.patient,
        "t1ce_date": d["t1ce"], "flair_date": d["flair"], "t1_date": d["t1"],
        "consistent": len(uniq) == 1,
    })
df = pd.DataFrame(rows)
df.to_csv(r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\check_study_date_consistency.csv",
          index=False, encoding="utf-8-sig")

print("总例数:", len(df))
print("日期一致:", df.consistent.sum(), " 不一致:", (~df.consistent).sum())
print()
print("=== 不一致病例 ===")
print(df[~df.consistent].to_string(index=False))
print()
print("=== 无日期层（结构不同）病例数 ===")
print((df.t1ce_date == "<no-date-layer>").sum())
print(df[df.t1ce_date == "<no-date-layer>"].patient.tolist())
