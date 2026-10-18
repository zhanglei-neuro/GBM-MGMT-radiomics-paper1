# -*- coding: utf-8 -*-
r"""
诊断：series_selection.csv（生产流程的序列选择）在 E 盘的可解析性。
"""
import pandas as pd, os

SEL = r"C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\series_selection.csv"
COHORT = r"D:\Raidionics_data\input\2"

s = pd.read_csv(SEL)
s = s[s.patient.notna()]
pats = set(s.patient.astype(str))
coh = set(os.listdir(COHORT))
print("队列 109 例中 series_selection 未覆盖:", sorted(coh - pats))
print("series_selection 中非队列病例:", sorted(pats - coh))
print()

ok = tot = 0
bad = []
for _, r in s.iterrows():
    for c in ["t1ce", "flair", "t1"]:
        p = str(r[c])
        if not p or p == "nan":
            bad.append((r.patient, c, "<EMPTY>"))
            continue
        p2 = p.replace("D:\\", "E:\\").replace("D:/", "E:/")
        tot += 1
        if os.path.exists(p2):
            ok += 1
        else:
            bad.append((r.patient, c, p2))
print(f"E盘可解析: {ok}/{tot}")
print("\n失败清单:")
for b in bad:
    print("  ", b)
