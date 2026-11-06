# -*- coding: utf-8 -*-
r"""Diagnostics: whether series_selection.csv (the series chosen by the production pipeline) resolves on the E: drive.
"""
import pandas as pd, os

SEL = r"C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\series_selection.csv"
COHORT = r"D:\Raidionics_data\input\2"

s = pd.read_csv(SEL)
s = s[s.patient.notna()]
pats = set(s.patient.astype(str))
coh = set(os.listdir(COHORT))
print("Cohort cases (of 109) not covered by series_selection:", sorted(coh - pats))
print("Cases in series_selection that are not in the cohort:", sorted(pats - coh))
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
print(f"Resolvable on the E: drive: {ok}/{tot}")
print("\nFailure list:")
for b in bad:
    print("  ", b)
