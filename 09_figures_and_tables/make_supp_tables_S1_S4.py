# -*- coding: utf-8 -*-
"""Build Supplementary Tables S1 (acquisition parameters) and S2 (radiomics settings)
for Paper 1.

Data sources (all primary, no re-parsing of DICOM geometry):
  * input/2/<patient>/T0/*.nii.gz     — production NIfTI volumes written by
    convert_inputs.py ("pure copy" DICOM->NIfTI, no resampling) => true acquisition
    geometry, read directly with nibabel.
  * series_selection.csv              — DICOM SeriesDescription + study date per role.
  * dicom_inventory.csv               — fallback for cases not in series_selection.
  * GBM_clinical_data.xlsx            — MRI_vendor.
  * retained source DICOM (8 cases)   — full header incl. TR/TE/field strength (S1c).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name   # noqa: E402  患者编号 <-> 本地姓名
import os, re, csv, glob
import numpy as np
import pandas as pd
import nibabel as nib
import pydicom

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side

BASE = r"C:/Users/zhanglei/WorkBuddy/2026-09-18-22-15-57/scripts"
SEL  = os.path.join(BASE, "series_selection.csv")
INV  = os.path.join(BASE, "dicom_inventory.csv")
XLS  = r"D:/1桌面文件/AI/GBM_clinical_data.xlsx"
INP  = r"D:/Raidionics_data/input/2"
OUTD = r"D:/WorkBuddy/2026-10-06-10-35-00/outputs/paper1/supplementary/tables"
BS   = chr(92)
os.makedirs(OUTD, exist_ok=True)

# 临床表 Name 列的拼写变体由 tools/patient_alias_local.csv 提供
# (本脚本不使用该表; build_acquisition_table.py 通过 patient_ids.resolve() 使用)

sel = {r["patient"]: r for r in csv.DictReader(open(SEL, encoding="utf-8-sig"))}
inv = list(csv.DictReader(open(INV, encoding="utf-8-sig")))
def nrm_p(p):
    return str(p).replace("/", BS).strip().rstrip(BS)


inv_map = {(r["patient"], nrm_p(r["series_dir"])): r for r in inv}
inv_by_pat = {}
for r in inv:
    inv_by_pat.setdefault(r["patient"], []).append(r)

# vendor: 1) vendor_map.csv (built from original DICOM series naming, most authoritative)
#         2) retained DICOM Manufacturer for add-on cases
vm = pd.read_csv(os.path.join(BASE, "vendor_map.csv"), encoding="utf-8-sig")
vendor_map = dict(zip(vm["patient_pinyin"].astype(str).str.strip(),
                      vm["MRI_vendor"].astype(str).str.strip()))
clin = pd.read_excel(XLS, sheet_name=0)   # kept for reference / cohort cross-check
clin["Name"] = clin["Name"].astype(str).str.strip()

# 加测批病例（DICOM 保留在 Raidionics_data/disc_* 下）; 键为本地目录名
ADDON = {to_name(x) for x in
         ["P001", "P004", "P025", "P057", "P060", "P067", "P089"]}


def parse_desc(s):
    m = re.match(r"(.+?)\s*\((\d+)slices\s+sp=([\d.]+)\)", str(s))
    if m:
        return m.group(1).strip(), int(m.group(2)), float(m.group(3))
    return str(s), None, None


def find_file(d, kw):
    if not os.path.isdir(d):
        return None
    for f in sorted(os.listdir(d)):
        if kw in f.lower():
            return os.path.join(d, f)
    return None


# DICOM SeriesDescription fallback for the add-on cases retained on disc
DISC_ROOT = r"D:/Raidionics_data"
DISC = {to_name(x): os.path.join(DISC_ROOT, "disc_" + to_name(x).lower())
        for x in ["P001", "P004", "P025", "P057", "P060", "P067", "P089"]}
disc_desc = {}
for pat, root in DISC.items():
    if not os.path.isdir(root):
        continue
    for ser in sorted(os.listdir(root)):
        sdir = os.path.join(root, ser)
        dcms = sorted(glob.glob(os.path.join(sdir, "*.dcm")))
        if not dcms:
            continue
        try:
            d = pydicom.dcmread(dcms[0], stop_before_pixels=True)
        except Exception:
            continue
        desc = str(getattr(d, "SeriesDescription", ""))
        # role inference
        dl = desc.lower()
        role = ("t1ce" if ("+" in desc or "+c" in dl or "gd" in dl)
                else "flair" if ("flair" in dl or "dark-fluid" in dl)
                else "t1" if dl.startswith("t1") else "other")
        if role != "other":
            disc_desc.setdefault(pat, {})[role] = desc

ROLES = [("T1c (contrast-enhanced T1)", "t1ce", "t1-gd"),
         ("T1w (non-contrast T1)", "t1", "t1w"),
         ("FLAIR", "flair", "t2-flair")]

rows = []
pats = sorted(p for p in os.listdir(INP)
              if os.path.isdir(os.path.join(INP, p, "T0")))
for pat in pats:
    vend = vendor_map.get(pat)
    batch = "main (batch 2)" if pat not in ADDON else "add-on (retained DICOM)"
    srow = sel.get(pat)
    for label, role, kw in ROLES:
        f = find_file(os.path.join(INP, pat, "T0"), kw)
        if f is None:
            continue
        im = nib.load(f)
        nx, ny, nz = im.shape
        sx, sy, sz = (float(x) for x in im.header.get_zooms()[:3])
        # description & study date
        desc = study = ""
        if srow is not None:
            desc, _, _ = parse_desc(srow.get(role + "_desc", ""))
            key = nrm_p(srow.get(role, ""))
            invr = inv_map.get((pat, key))
            if invr:
                study = invr["study_date"]
        if not desc and pat in disc_desc and role in disc_desc[pat]:
            desc = disc_desc[pat][role]
        if not study:
            for r in inv_by_pat.get(pat, []):
                study = r["study_date"]
                break
        rows.append({
            "Patient": pat,
            "Cohort batch": batch,
            "MRI vendor": vend if vend else
                ("Siemens" if pat in ADDON else "n/a"),
            "Sequence role": label,
            "DICOM SeriesDescription": desc,
            "Study date": study,
            "Matrix (rows x cols)": "%d x %d" % (nx, ny),
            "Rows": nx, "Cols": ny, "Slices": nz,
            "In-plane pixel spacing (mm)": round(sx, 4),
            "Slice spacing (mm)": round(sz, 4),
            "FOV x (mm)": round(nx * sx, 1),
            "FOV y (mm)": round(ny * sy, 1),
            "Contrast agent": "Gd-DTPA" if role == "t1ce" else "",
        })

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUTD, "S1b_acquisition_per_case.csv"),
          index=False, encoding="utf-8-sig")


# ---------------------------------------------------------------- S1a summary
def rng(s, fmt="%.2f"):
    s = s.dropna()
    if len(s) == 0:
        return ""
    return (fmt + " [" + fmt + "-" + fmt + "]") % (s.median(), s.min(), s.max())


def matrix_rng(g):
    a = g["Matrix (rows x cols)"].value_counts()
    if len(a) == 1:
        return a.index[0]
    return "; ".join("%s (n=%d)" % (k, v) for k, v in a.head(3).items())


def prot_rng(g):
    a = g["DICOM SeriesDescription"].replace("", np.nan).dropna().value_counts()
    if len(a) == 0:
        return ""
    return "; ".join("%s (n=%d)" % (k, v) for k, v in a.head(3).items())


sum_rows = []
for vend in ["Siemens", "GE", "UIH"]:
    for label, role, kw in ROLES:
        g = df[(df["MRI vendor"] == vend) & (df["Sequence role"] == label)]
        if len(g) == 0:
            continue
        sum_rows.append({
            "Vendor": vend, "Sequence": label, "n patients": g["Patient"].nunique(),
            "Acquisition protocol(s)": prot_rng(g),
            "Matrix": matrix_rng(g),
            "In-plane pixel spacing (mm)": rng(g["In-plane pixel spacing (mm)"]),
            "Slice spacing (mm)": rng(g["Slice spacing (mm)"]),
            "Slices": rng(g["Slices"].astype(float), fmt="%.0f"),
            "FOV x (mm)": rng(g["FOV x (mm)"], fmt="%.0f"),
        })
s1a = pd.DataFrame(sum_rows)
s1a.to_csv(os.path.join(OUTD, "S1a_acquisition_by_vendor.csv"),
           index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- S1c retained DICOM
ROOTS = list(DISC.values()) + [os.path.join(r"D:/NeuroImagepre", to_name("P001"))]
s1c_rows = []
for root in ROOTS:
    if not os.path.isdir(root):
        continue
    seen = set()
    for dp, dns, fns in os.walk(root):
        dcms = sorted(x for x in fns if x.lower().endswith(".dcm"))
        if not dcms:
            continue
        try:
            d = pydicom.dcmread(os.path.join(dp, dcms[0]), stop_before_pixels=True)
        except Exception:
            continue
        pat = os.path.basename(root).replace("disc_", "")
        ser = os.path.basename(dp)
        if (pat, ser) in seen:
            continue
        seen.add((pat, ser))
        ps = d.PixelSpacing if hasattr(d, "PixelSpacing") else [None, None]
        s1c_rows.append({
            "Patient": pat, "Series dir": ser,
            "Manufacturer": str(getattr(d, "Manufacturer", "")),
            "Model": str(getattr(d, "ManufacturerModelName", "")),
            "Field strength (T)": getattr(d, "MagneticFieldStrength", ""),
            "SeriesDescription": str(getattr(d, "SeriesDescription", "")),
            "TR (ms)": getattr(d, "RepetitionTime", ""),
            "TE (ms)": getattr(d, "EchoTime", ""),
            "TI (ms)": getattr(d, "InversionTime", ""),
            "Flip angle (deg)": getattr(d, "FlipAngle", ""),
            "Slice thickness (mm)": getattr(d, "SliceThickness", ""),
            "Slice spacing (mm)": getattr(d, "SpacingBetweenSlices", ""),
            "Pixel spacing (mm)": round(float(ps[0]), 4) if ps and ps[0] else None,
            "Matrix": "%d x %d" % (getattr(d, "Rows", 0), getattr(d, "Cols", 0)),
            "Acquisition type": str(getattr(d, "MRAcquisitionType", "")),
            "Software": str(getattr(d, "SoftwareVersions", "")),
            "Study date": str(getattr(d, "StudyDate", "")),
        })
s1c = pd.DataFrame(s1c_rows)
s1c.to_csv(os.path.join(OUTD, "S1c_retained_DICOM_parameters.csv"),
           index=False, encoding="utf-8-sig")

# ---------------------------------------------------------------- S2 settings
s2 = pd.DataFrame([
    ("Software", "PyRadiomics 3.1.0 (Python 3.9, Windows; dedicated conda env)"),
    ("Extraction script", "extract_features.py — deterministic production script, "
                          "verbatim source provided with the supplement"),
    ("Image types", "Original; Wavelet (8 decompositions: LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH)"),
    ("Feature classes", "shape (3D, primary pairs only); firstorder (18); glcm (24); glrlm (16); "
                        "glszm (16); gldm (14); ngtdm (5)"),
    ("binWidth", "25 (fixed bin width)"),
    ("Resampling", "None — features computed on the native voxel grid of each sequence"),
    ("Intensity normalisation", "None — original DICOM-derived intensities retained "
                                "(inter-vendor harmonisation evaluated separately in the main text)"),
    ("Interpolator", "n/a (no resampling during extraction)"),
    ("Label value", "1"),
    ("correctMask", "True (residual geometry mismatches are auto-resampled rather than dropping the case)"),
    ("geometryTolerance", "1e-3"),
    ("force2D", "False (full 3D texture matrices)"),
    ("Neighbourhood / distance", "PyRadiomics defaults (offset distance = 1 voxel)"),
    ("ROI x sequence pairing (primary)",
     "TumorCore x T1c; Necrosis x T1c; Edema x FLAIR; FLAIRChanges x FLAIR — "
     "Original + Wavelet + shape: 93 x 9 + 14 = 851 features per pair; 4 pairs = 3404"),
    ("ROI x sequence pairing (cross)",
     "Remaining 8 ROI x sequence combinations — Original only: 93 features each = 744"),
    ("Total features per patient", "4148 (3404 primary + 744 cross)"),
    ("Segmentation source",
     "Raidionics (v1.2) neuro pre-operative pipeline: MRI_TumorCore on T1c, dedicated t1c-based "
     "necrosis model, MRI_FLAIRChanges on FLAIR, plus brain masking; no manual edits"),
    ("Multi-sequence alignment",
     "T1c is the reference grid. T1w/FLAIR with native geometry equal to T1c are used as-is; "
     "otherwise resampled with the Raidionics ANTs forward transform (affine + warp, linear "
     "interpolation), cached per patient. Masks are matched to the per-pair image grid "
     "(nearest-neighbour resample fallback)."),
    ("Necrosis mask note",
     "Necrosis masks were produced by a second dedicated Raidionics run (necrosis_t1c_config.ini, "
     "custom model folder models_nec_t1c) after the initial batch, to correct known necrosis "
     "under-segmentation; these masks override the originals in patients_nec2."),
], columns=["Item", "Value"])

# ---------------------------------------------------------------- workbook
XL_H = "FF4472C4"
thin = Side(style="thin", color="FFBFBFBF")
bd = Border(left=thin, right=thin, top=thin, bottom=thin)


def style(ws, ncol, nrow, widths, wrap_cols=()):
    for j in range(1, ncol + 1):
        c = ws.cell(row=1, column=j)
        c.font = Font(bold=True, color="FFFFFFFF")
        c.fill = PatternFill("solid", fgColor=XL_H)
        c.alignment = Alignment(horizontal="center", vertical="center", wrapText=True)
        c.border = bd
    for i in range(2, nrow + 1):
        for j in range(1, ncol + 1):
            cc = ws.cell(row=i, column=j)
            cc.border = bd
            cc.alignment = Alignment(vertical="center", wrapText=(j in wrap_cols))
    for col, w in widths.items():
        ws.column_dimensions[col].width = w
    ws.freeze_panes = "A2"


def dump(wb, sheet, data, widths, wrap=()):
    ws = wb.create_sheet(sheet)
    ws.append(list(data.columns))
    for _, r in data.iterrows():
        ws.append([None if (isinstance(v, float) and np.isnan(v)) else v
                   for v in r.tolist()])
    style(ws, len(data.columns), len(data) + 1, widths, wrap_cols=wrap)
    return ws


wb = Workbook()
ws = wb.active
ws.title = "README"
readme = [
    ["Supplementary tables — Paper 1 (multi-vendor Raidionics radiomics, MGMT prediction)"],
    [""],
    ["Table", "Content", "Rows"],
    ["S1a", "MRI acquisition parameters summarised by vendor and sequence role", len(s1a)],
    ["S1b", "Per-case acquisition parameters (all sequences entering the pipeline)", len(df)],
    ["S1c", "Full DICOM header for the retained source DICOM subset (n=%d cases, all Siemens 3T)" % s1c["Patient"].nunique(), len(s1c)],
    ["S2", "PyRadiomics feature extraction settings (production configuration)", len(s2)],
    [""],
    ["DATA PROVENANCE"],
    ["Geometry", "Read directly from the production NIfTI volumes (input/2), which were converted "
                 "from the original pre-operative DICOM with a pure-copy converter (no resampling, "
                 "no skew correction) — therefore identical to the acquisition geometry."],
    ["Protocol names / dates", "series_selection.csv and dicom_inventory.csv, built when the original "
                                "DICOM tree was still on-line."],
    ["Vendor labels", "GBM_clinical_data.xlsx column MRI_vendor; add-on cases labelled from their "
                      "retained DICOM headers."],
    ["Radiomics settings", "scripts/extract_features.py — the exact production script that generated "
                           "radiomics_features_long.csv (4148 features x 106 patients)."],
    [""],
    ["KNOWN LIMITATION — field strength / TR / TE"],
    ["The pre-operative source DICOM tree for the main cohort was archived off-line after conversion; "
     "only %d cases retain local DICOM (all Siemens 3T: Verio, MAGNETOM Vida)." % s1c["Patient"].nunique()],
    ["Consequently MagneticFieldStrength / RepetitionTime / EchoTime cannot be reported for the full "
     "cohort from DICOM. Table S1c documents what the retained subset provides; if the archived "
     "pre-operative DICOM (or a PACS export) becomes available, S1b can be extended with these "
     "columns in a single pass."],
    ["Reference only: post-operative DICOM (n=64 patients, same centre, overlapping period) shows "
     "GE SIGNA Premier (n=26), Siemens Verio (n=31), Siemens MAGNETOM Vida (n=5) and UIH uMR560 "
     "(n=2); field strength 3T in 62/64 and 1.5T in 2/64 (both UIH)."],
    [""],
    ["Generated", pd.Timestamp.now().strftime("%Y-%m-%d %H:%M")],
]
for row in readme:
    ws.append(row)
ws.column_dimensions["A"].width = 28
ws.column_dimensions["B"].width = 130
ws.column_dimensions["C"].width = 10
ws["A1"].font = Font(bold=True, size=12)

dump(wb, "S1a_by_vendor", s1a,
     {"A": 10, "B": 26, "C": 10, "D": 56, "E": 30, "F": 16, "G": 13, "H": 14, "I": 12},
     wrap=(4, 5))
dump(wb, "S1b_per_case", df,
     {"A": 16, "B": 20, "C": 10, "D": 24, "E": 26, "F": 11, "G": 14,
      "H": 7, "I": 7, "J": 8, "K": 14, "L": 12, "M": 10, "N": 10, "O": 12})
dump(wb, "S1c_retained_DICOM", s1c,
     {"A": 14, "B": 18, "C": 12, "D": 16, "E": 10, "F": 24, "G": 9, "H": 9, "I": 9,
      "J": 10, "K": 11, "L": 11, "M": 12, "N": 11, "O": 12, "P": 14, "Q": 11})
dump(wb, "S2_radiomics_settings", s2, {"A": 32, "B": 125}, wrap=(2,))

out_xlsx = os.path.join(OUTD, "Supplementary_Tables_S1_S2.xlsx")
wb.save(out_xlsx)
print("saved:", out_xlsx)
print("S1a rows:", len(s1a), " S1b rows:", len(df), " S1c rows:", len(s1c))
print("\n=== S1a preview ===")
print(s1a.to_string(index=False, max_colwidth=44))
print("\n=== vendor x batch counts (per patient) ===")
print(df.groupby(["MRI vendor", "Cohort batch"])["Patient"].nunique())
