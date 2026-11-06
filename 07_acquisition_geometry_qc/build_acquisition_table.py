# -*- coding: utf-8 -*-
r"""
Whole-database MRI acquisition audit -- NOT a supplementary table
Data sources: E:\NeuroImagepre (107 cases, including 3 with a PAT/STUDY/SERIES structure)
              + D:\Raidionics_data\disc_* (2 cases: P004 / P025)

This script works on the full 109-case database and all three series of every case, so its output
includes the single United Imaging 1.5 T case and does not match the 89-case paper cohort. It is
kept for traceability of the vendor/geometry audit only. The tables actually submitted as
Supplementary Tables S1a-S1c are produced by `build_s1_paper1_used.py`, which keeps one series per
sequence for the 89 cases; the files written here were withdrawn from the submission package on
2026-11-06 (see supplementary_tables/TABLE_LEGENDS.md).

Output:
  S1a_acquisition_by_vendor.csv   summary by vendor x sequence role
  S1b_acquisition_per_case.csv    per case per sequence detail (with TR/TE/TI/field strength/model)
  S1c_vendor_signature.csv        vendor technical signature (field strength and model distribution)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name, resolve   # noqa: E402  patient code <-> local name
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

# ---------- the 6 cases not covered by series_selection: specified by hand (from ContrastBolusAgent and the series name) ----------
# value = (root type, study_dir, series_dir)
# keys are local patient directory names (recovered from the de-identified code via to_name())
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
    """Extract the number; returns None on failure"""
    if v is None or v == "":
        return None
    s = str(v).replace("\\", " ")
    try:
        out = [float(x) for x in s.split() if x not in ("", "nan")]
        return out[0] if out else None
    except Exception:
        return None


def valstr(v):
    """Convert a DICOM value to a string consistently (MultiValue -> space separated, no quotes)"""
    if v is None:
        return ""
    if isinstance(v, (MultiValue, list, tuple)):
        return " ".join(valstr(x) for x in v)
    return str(v).replace("\\", " ").strip()


def read_header(series_dir):
    """Read the header of the first DICOM in the series directory"""
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
    """Locate series_dir under root (restricted to the study level)"""
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
            # formatting
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

# ---------- QC: compare with the date of surgery ----------
    xls = r"D:\1桌面文件\AI\GBM_clinical_data.xlsx"
    x = pd.read_excel(xls)
    x.columns = [str(c).strip() for c in x.columns]
# Name column of the clinical table -> local patient directory name: first through the spelling-variant table and code normalisation, then back to a directory name
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

# whether the three sequences of each case share one date
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
    print("\n=== QC flag summary (by case) ===")
    q = df.groupby("Patient")["QC flag"].first().value_counts()
    print(q.to_string())
    print("\n=== cases that are not OK ===")
    badp = df.groupby("Patient")["QC flag"].first()
    badp = badp[badp != "OK"]
    for p, f in badp.items():
        sub = df[df.Patient == p]
        print(f"  {p:14s} {f:18s} " + "; ".join(
            f"{r['Sequence role']}={r.StudyDate}(Δ{r['Days from surgery']:+.0f}d)" for _, r in sub.iterrows()))

    if unmatched:
        print("\nunmatched list:")
        for u in unmatched:
            print("  ", u)

# ---------- S1a: summary by vendor x sequence role ----------
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

# ---------- S1c: vendor technical signature (field strength x model) ----------
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

# save the full detail for checking
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
