# -*- coding: utf-8 -*-
r"""
Paper 1 补充材料 Table S1 — 全队列 MRI 采集参数重建
数据源：E:\NeuroImagepre（107 例，含 3 例 PAT/STUDY/SERIES 结构）
        + D:\Raidionics_data\disc_*（2 例：P004 / P025）

输出：
  S1a_acquisition_by_vendor.csv   按厂商 × 序列角色汇总
  S1b_acquisition_per_case.csv    逐例逐序列明细（含 TR/TE/TI/场强/机型）
  S1c_vendor_signature.csv        厂商技术参数签名（场强/机型分布）
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name, resolve   # noqa: E402  患者编号 <-> 本地姓名
import os, csv, json, re
import numpy as np
import pandas as pd
import pydicom
from pydicom.multival import MultiValue

OUT = r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\tables"
SCR = r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts"
os.makedirs(OUT, exist_ok=True)

SEL = r"C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\series_selection.csv"
COHORT_DIR = r"D:\Raidionics_data\input\2"
E_ROOT = r"E:\NeuroImagepre"
DISC = r"D:\Raidionics_data"

ROLES = ["T1c", "T1w", "FLAIR"]

# ---------- 6 例未被 series_selection 覆盖的病例：手工指定（依据 ContrastBolusAgent 与序列名） ----------
# 值 = (根目录类型, study_dir, series_dir)
# 键为本地患者目录名（由去标识编号经 to_name() 还原）
MANUAL = {
    to_name("P001"): {"study": "STD00001",    "root": "E",
                      "T1c": "SER00008",       "T1w": "SER00004",       "FLAIR": "SER00007"},
    to_name("P057"): {"study": "PAT00000001", "root": "E",
                      "T1c": "SERIES00000012",  "T1w": "SERIES00000004", "FLAIR": "SERIES00000007"},
    to_name("P060"): {"study": "STD00002",    "root": "E",
                      "T1c": "SER00008",       "T1w": "SER00004",       "FLAIR": "SER00007"},
    to_name("P089"): {"study": "PAT00000001", "root": "E",
                      "T1c": "SERIES00000009",  "T1w": "SERIES00000004", "FLAIR": "SERIES00000005"},
    to_name("P004"): {"study": "disc",        "root": "D",
                      "T1c": "SERIES00000008",  "T1w": "SERIES00000004", "FLAIR": "SERIES00000007"},
    to_name("P025"): {"study": "disc",        "root": "D",
                      "T1c": "SER00011",       "T1w": "SER00007",       "FLAIR": "SER00008"},
}

HDR_FIELDS = [
    ("Manufacturer", "Manufacturer"),
    ("ManufacturerModelName", "Model"),
    ("MagneticFieldStrength", "FieldStrength_T"),
    ("SoftwareVersions", "Software"),
    ("SeriesDescription", "SeriesDescription"),
    ("RepetitionTime", "TR_ms"),
    ("EchoTime", "TE_ms"),
    ("InversionTime", "TI_ms"),
    ("FlipAngle", "FlipAngle_deg"),
    ("MRAcquisitionType", "AcqType"),
    ("SliceThickness", "SliceThickness_mm"),
    ("SpacingBetweenSlices", "SliceSpacing_mm"),
    ("Rows", "Rows"),
    ("Columns", "Cols"),
    ("PixelSpacing", "PixelSpacing_mm"),
    ("ContrastBolusAgent", "ContrastAgent"),
    ("StudyDate", "StudyDate"),
    ("SeriesNumber", "SeriesNumber"),
]


def num(v):
    """提取数值，失败返回 None"""
    if v is None or v == "":
        return None
    s = str(v).replace("\\", " ")
    try:
        out = [float(x) for x in s.split() if x not in ("", "nan")]
        return out[0] if out else None
    except Exception:
        return None


def valstr(v):
    """统一把 DICOM 值转字符串（MultiValue -> 空格分隔，不带引号）"""
    if v is None:
        return ""
    if isinstance(v, (MultiValue, list, tuple)):
        return " ".join(valstr(x) for x in v)
    return str(v).replace("\\", " ").strip()


def read_header(series_dir):
    """读序列目录第一个 DICOM 的头"""
    if not os.path.isdir(series_dir):
        return None
    files = sorted(os.listdir(series_dir))
    for fn in files:
        fp = os.path.join(series_dir, fn)
        if not os.path.isfile(fp):
            continue
        try:
            ds = pydicom.dcmread(fp, stop_before_pixels=True, force=True)
        except Exception:
            continue
        rec = {}
        for tag, name in HDR_FIELDS:
            rec[name] = valstr(getattr(ds, tag, None))
        rec["n_slices"] = len(files)
        return rec
    return None


def find_series_dir(root, study, series_name):
    """在 root 下定位 series_dir（限定 study 层）"""
    for dirpath, dirnames, filenames in os.walk(root):
        if os.path.basename(dirpath) != series_name:
            continue
        if study and study != "disc" and study not in dirpath:
            continue
        if filenames:
            return dirpath
    return None


def main():
    allc = sorted(os.listdir(COHORT_DIR))
    sel = pd.read_csv(SEL)
    sel = sel[sel.patient.notna()].set_index("patient")

    rows_s1b = []
    unmatched = []

    for pat in allc:
        recs = {}
        if pat in MANUAL:
            m = MANUAL[pat]
            root = os.path.join(E_ROOT, pat) if m["root"] == "E" else os.path.join(DISC, "disc_" + pat.lower())
            for role in ROLES:
                sd = find_series_dir(root, m["study"], m[role])
                if sd is None:
                    unmatched.append((pat, role, "dir-not-found"))
                    continue
                h = read_header(sd)
                if h is None:
                    unmatched.append((pat, role, "no-dicom"))
                    continue
                recs[role] = (sd, h)
        elif pat in sel.index:
            r = sel.loc[pat]
            colmap = {"T1c": "t1ce", "T1w": "t1", "FLAIR": "flair"}
            for role in ROLES:
                p = str(r[colmap[role]])
                if not p or p == "nan":
                    unmatched.append((pat, role, "empty-in-selection"))
                    continue
                p2 = p.replace("D:\\", "E:\\").replace("D:/", "E:/")
                h = read_header(p2)
                if h is None:
                    unmatched.append((pat, role, "not-readable:" + p2))
                    continue
                recs[role] = (p2, h)
        else:
            unmatched.append((pat, "-", "patient-not-in-selection"))
            continue

        for role in ROLES:
            if role not in recs:
                continue
            sd, h = recs[role]
            row = {
                "Patient": pat,
                "DICOM source": "E:\\NeuroImagepre" if sd.upper().startswith("E:") else "D:\\Raidionics_data",
                "Sequence role": role,
                "Series dir": os.path.relpath(sd, E_ROOT if sd.upper().startswith("E:") else DISC),
            }
            row.update(h)
            # 格式化
            for k in ["TR_ms", "TE_ms", "TI_ms", "FlipAngle_deg", "SliceThickness_mm", "SliceSpacing_mm"]:
                v = num(h.get(k))
                row[k] = "" if v is None else (f"{v:g}")
            r_, c_ = num(h.get("Rows")), num(h.get("Cols"))
            row["Matrix"] = f"{int(r_)} x {int(c_)}" if r_ and c_ else ""
            ps = h.get("PixelSpacing_mm", "")
            ps_parts = re.findall(r"\d+\.?\d*", str(ps))
            row["Pixel spacing row (mm)"] = f"{float(ps_parts[0]):.4g}" if len(ps_parts) >= 1 else ""
            row["Pixel spacing col (mm)"] = f"{float(ps_parts[1]):.4g}" if len(ps_parts) >= 2 else ""
            row["FOV x (mm)"] = f"{c_ * float(ps_parts[1]):.1f}" if (c_ and len(ps_parts) >= 2) else ""
            row["FOV y (mm)"] = f"{r_ * float(ps_parts[0]):.1f}" if (r_ and len(ps_parts) >= 1) else ""
            fs = num(h.get("FieldStrength_T"))
            row["FieldStrength_T"] = "" if fs is None else f"{fs:g}"
            vendor = str(h.get("Manufacturer", "")).strip()
            row["Vendor"] = normalize_vendor(vendor)
            rows_s1b.append(row)

    df = pd.DataFrame(rows_s1b)

    # ---------- QC：与手术日期比对 ----------
    xls = r"D:\1桌面文件\AI\GBM_clinical_data.xlsx"
    x = pd.read_excel(xls)
    x.columns = [str(c).strip() for c in x.columns]
    # 临床表 Name 列 → 本地患者目录名：先经拼写变体表与编号归一, 再还原为目录名
    x["pid"] = (x["Name"].astype(str).str.replace(r"\(.*?\)", "", regex=True).str.strip()
                .map(lambda n: to_name(resolve(n))))
    surg = dict(zip(x.pid, pd.to_datetime(x["Surgery date"], errors="coerce")))

    def to_dt(v):
        try:
            return pd.Timestamp(str(int(float(v))))
        except Exception:
            return pd.NaT

    df["StudyDate_dt"] = df.StudyDate.map(to_dt)
    df["SurgeryDate"] = df.Patient.map(surg)
    df["Days from surgery"] = (df.StudyDate_dt - df.SurgeryDate).dt.days

    # 每例三序列日期是否一致
    cons = df.groupby("Patient").StudyDate.nunique().rename("n_unique_dates")
    df = df.merge(cons, on="Patient", how="left")

    def qc(row):
        d = row["Days from surgery"]
        if pd.isna(d):
            return "NO_SURGERY_DATE"
        if row.n_unique_dates > 1:
            return "MIXED_STUDY_DATES"
        if d > 0:
            return "POSTOP_IMAGE"
        return "OK"

    df["QC flag"] = df.apply(qc, axis=1)
    df = df.drop(columns=["StudyDate_dt", "SurgeryDate"])

    df.to_csv(os.path.join(OUT, "S1b_acquisition_per_case.csv"), index=False, encoding="utf-8-sig")
    print("S1b rows:", len(df), "patients:", df.Patient.nunique())
    print("\n=== QC flag 汇总（按例）===")
    q = df.groupby("Patient")["QC flag"].first().value_counts()
    print(q.to_string())
    print("\n=== 非 OK 病例 ===")
    badp = df.groupby("Patient")["QC flag"].first()
    badp = badp[badp != "OK"]
    for p, f in badp.items():
        sub = df[df.Patient == p]
        print(f"  {p:14s} {f:18s} " + "; ".join(
            f"{r['Sequence role']}={r.StudyDate}(Δ{r['Days from surgery']:+.0f}d)" for _, r in sub.iterrows()))

    if unmatched:
        print("\n未匹配清单:")
        for u in unmatched:
            print("  ", u)

    # ---------- S1a: 按厂商 × 序列角色汇总 ----------
    def rng(s):
        v = pd.to_numeric(s, errors="coerce").dropna()
        if v.empty:
            return ""
        return f"{v.median():g} [{v.min():g}-{v.max():g}]"

    s1a = []
    for (ven, role), g in df.groupby(["Vendor", "Sequence role"]):
        prot = g.SeriesDescription.value_counts()
        prot_s = "; ".join(f"{k} (n={v})" for k, v in prot.items())
        mat = g.Matrix.value_counts()
        mat_s = "; ".join(f"{k} (n={v})" for k, v in mat.items())
        fs = g.FieldStrength_T.value_counts()
        fs_s = "; ".join(f"{k}T (n={v})" for k, v in fs.items())
        mdl = g.Model.value_counts()
        mdl_s = "; ".join(f"{k} (n={v})" for k, v in mdl.items())
        s1a.append({
            "Vendor": ven, "Sequence": role, "n patients": g.Patient.nunique(),
            "Field strength": fs_s, "Scanner model": mdl_s,
            "Acquisition protocol(s)": prot_s,
            "TR (ms)": rng(g.TR_ms), "TE (ms)": rng(g.TE_ms), "TI (ms)": rng(g.TI_ms),
            "Flip angle (deg)": rng(g.FlipAngle_deg),
            "Matrix": mat_s,
            "In-plane pixel spacing (mm)": rng(g["Pixel spacing col (mm)"]),
            "Slice thickness (mm)": rng(g.SliceThickness_mm),
            "Slice spacing (mm)": rng(g.SliceSpacing_mm),
            "Slices": rng(g.n_slices),
        })
    a = pd.DataFrame(s1a).sort_values(["Vendor", "Sequence"])
    a.to_csv(os.path.join(OUT, "S1a_acquisition_by_vendor.csv"), index=False, encoding="utf-8-sig")
    print("\nS1a rows:", len(a))
    print(a[["Vendor", "Sequence", "n patients", "Field strength", "TR (ms)", "TE (ms)", "TI (ms)"]].to_string(index=False))

    # ---------- S1c: 厂商技术签名（场强 × 机型） ----------
    uniq = df.drop_duplicates(["Patient"])
    tech = uniq.groupby(["Vendor", "Model", "FieldStrength_T"]).agg(
        n_patients=("Patient", "nunique"),
    ).reset_index()
    soft = (uniq.groupby(["Vendor", "Model", "FieldStrength_T"])["Software"]
            .apply(lambda x: "; ".join(sorted({v for v in x if v and v != "nan"}))))
    tech = tech.merge(soft.rename("Software versions").reset_index(),
                      on=["Vendor", "Model", "FieldStrength_T"], how="left")
    tech.to_csv(os.path.join(OUT, "S1c_vendor_signature.csv"), index=False, encoding="utf-8-sig")
    print("\nS1c:")
    print(tech.to_string(index=False))

    # 保存完整明细供核查
    df.to_csv(os.path.join(SCR, "S1b_full_detailed.csv"), index=False, encoding="utf-8-sig")
    print("\nDONE")


def normalize_vendor(v):
    v = v.upper()
    if "SIEMENS" in v:
        return "Siemens"
    if "GE" in v:
        return "GE"
    if "UIH" in v or "UNITED IMAGING" in v:
        return "UIH"
    if "PHILIPS" in v:
        return "Philips"
    return v.strip() or "Unknown"


if __name__ == "__main__":
    main()
