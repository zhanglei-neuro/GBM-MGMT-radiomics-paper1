# -*- coding: utf-8 -*-
r"""
补充扫描：E:\NeuroImagepre 中 3 个 PAT/STUDY/SERIES/IMG 结构的病例
+ D:\Raidionics_data\disc_* 的 2 例（P004 / P025）。
统一逻辑：找出所有"直接包含 DICOM 文件的目录"，每个目录视为一个 series。
追加到 pre_dicom_series_inventory.csv
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name   # noqa: E402  患者编号 <-> 本地姓名
import os, csv
import pydicom
import pandas as pd

TAGS = [
    "SeriesDescription", "StudyDescription", "ProtocolName",
    "Manufacturer", "ManufacturerModelName",
    "MagneticFieldStrength", "FieldStrength",
    "TR", "RepetitionTime", "TE", "EchoTime", "InversionTime", "FlipAngle",
    "SliceThickness", "SpacingBetweenSlices", "PixelSpacing",
    "Rows", "Columns", "AcquisitionMatrix", "AcquisitionType",
    "ContrastBolusAgent", "ScanningSequence", "SequenceVariant",
    "MRAcquisitionType", "SoftwareVersions", "StudyDate", "SeriesNumber",
    "SliceLocation", "ImageOrientationPatient", "PatientName",
]

def is_dicom(fp):
    try:
        with open(fp, "rb") as f:
            f.seek(128)
            return f.read(4) == b"DICM"
    except Exception:
        return False

def probe(series_dir):
    files = [f for f in os.listdir(series_dir) if os.path.isfile(os.path.join(series_dir, f))]
    for fn in files:
        fp = os.path.join(series_dir, fn)
        if not is_dicom(fp):
            continue
        try:
            ds = pydicom.dcmread(fp, stop_before_pixels=True, force=True)
        except Exception:
            continue
        rec = {}
        for t in TAGS:
            v = getattr(ds, t, None)
            if v is None:
                rec[t] = ""
            elif isinstance(v, (list, tuple)):
                rec[t] = "\\".join(str(x) for x in v)
            else:
                rec[t] = str(v)
        rec["n_files"] = len(files)
        return rec
    return None

def walk_series(root):
    """返回 root 下所有直接包含 DICOM 的目录"""
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        if not filenames:
            continue
        f0 = os.path.join(dirpath, filenames[0])
        if is_dicom(f0):
            found.append(dirpath)
            dirnames[:] = []  # 不下钻
    return found

def main():
    rows = []
    targets = []
    # 1) E 盘三个特殊结构病例
    for p in [to_name(x) for x in ["P057", "P060", "P089"]]:
        targets.append((p, os.path.join(r"E:\NeuroImagepre", p)))
    # 2) disc_ 两例
    for p in [to_name(x) for x in ["P004", "P025"]]:
        targets.append((p, os.path.join(r"D:\Raidionics_data", "disc_" + p.lower())))

    for pat, root in targets:
        if not os.path.isdir(root):
            print(f"[SKIP] {pat}: {root} not found")
            continue
        sers = walk_series(root)
        print(f"[{pat}] {len(sers)} series dirs @ {root}")
        for sd in sers:
            rec = probe(sd)
            if rec is None:
                continue
            rec["patient"] = pat
            rec["study_date_dir"] = os.path.relpath(sd, root).split(os.sep)[0]
            rec["level1"] = os.path.basename(os.path.dirname(sd))
            rec["series_dir"] = os.path.basename(sd)
            rec["src_root"] = root
            rows.append(rec)

    out = r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\pre_dicom_series_inventory_extra.csv"
    cols = ["patient", "study_date_dir", "level1", "series_dir", "src_root"] + TAGS + ["n_files"]
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nDONE: {len(rows)} series -> {out}")

if __name__ == "__main__":
    main()
