# -*- coding: utf-8 -*-
r"""
全局核查：每例三个序列的 StudyDate vs 手术日期。
数据源：S1b_full_detailed.csv（含 DICOM StudyDate）+ GBM_clinical_data.xlsx（Surgery date）
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name, resolve   # noqa: E402  患者编号 <-> 本地姓名
import pandas as pd, numpy as np

DET = r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\S1b_full_detailed.csv"
XLS = r"D:\1桌面文件\AI\GBM_clinical_data.xlsx"

d = pd.read_csv(DET)
x = pd.read_excel(XLS)
x.columns = [str(c).strip() for c in x.columns]
x["Name"] = x["Name"].astype(str).str.strip()

# 别名（临床表 -> 队列拼音）
# 临床表 Name 列 → 本地患者目录名：先经拼写变体表与编号归一, 再还原为目录名
x["pid"] = (x["Name"].str.replace(r"\(.*?\)", "", regex=True).str.strip()
            .map(lambda n: to_name(resolve(n))))

surg = dict(zip(x.pid, pd.to_datetime(x["Surgery date"], errors="coerce")))
print("临床表可映射手术日期:", sum(v is not None and not pd.isna(v) for v in surg.values()), "/", len(x))

def to_dt(v):
    try:
        return pd.Timestamp(str(int(float(v))))
    except Exception:
        return pd.NaT

rows = []
for (pat, role), g in d.groupby(["Patient", "Sequence role"]):
    sd = to_dt(g.StudyDate.iloc[0]) if str(g.StudyDate.iloc[0]).strip() not in ("", "nan") else pd.NaT
    sur = surg.get(pat, pd.NaT)
    delta = (sd - sur).days if (pd.notna(sd) and pd.notna(sur)) else np.nan
    rows.append({"Patient": pat, "Role": role, "StudyDate": sd.date() if pd.notna(sd) else None,
                 "SurgeryDate": sur.date() if pd.notna(sur) else None,
                 "Delta_days": delta, "Postop?": (delta > 0) if pd.notna(delta) else None})

r = pd.DataFrame(rows)
r.to_csv(r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\check_study_vs_surgery.csv",
         index=False, encoding="utf-8-sig")

piv = r.pivot_table(index="Patient", columns="Role", values="Delta_days", aggfunc="first")
print("\n=== 序列日期晚于手术日期的记录（Postop）===")
bad = r[r["Postop?"] == True]
print(bad.to_string(index=False))
print(f"\n共 {bad.Patient.nunique()} 例受影响: {sorted(bad.Patient.unique())}")

print("\n=== 完整 Delta 分布（天）===")
print(piv.describe().to_string())
print("\n=== Delta > 0 的患者（任一序列）===")
print(piv[(piv > 0).any(axis=1)].to_string())
print("\n=== Delta 异常大（<-90 天，可能非术前最近一次）===")
print(piv[(piv < -90).any(axis=1)].to_string())
