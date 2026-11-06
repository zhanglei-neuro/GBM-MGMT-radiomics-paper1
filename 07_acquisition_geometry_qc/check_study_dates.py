# -*- coding: utf-8 -*-
r"""Batch check whether the three series in the production series selection (series_selection.csv) come from the same study.
Paths look like D:\NeuroImagepre\<pat>\<YYYYMMDD>\... or D:\NeuroImagepre\<pat>\ImageFileName_xxx\<n>
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

print("Total cases:", len(df))
print("Consistent study date:", df.consistent.sum(), " inconsistent:", (~df.consistent).sum())
print()
print("=== Inconsistent cases ===")
print(df[~df.consistent].to_string(index=False))
print()
print("=== Cases with no date level (different folder structure) ===")
print((df.t1ce_date == "<no-date-layer>").sum())
print(df[df.t1ce_date == "<no-date-layer>"].patient.tolist())
