# -*- coding: utf-8 -*-
r"""Scan the full preoperative DICOM archive under E:\NeuroImagepre and extract the acquisition parameters of every series.
Output: pre_dicom_series_inventory.csv (one row per patient per series)
"""
import os, sys, csv, glob, json
import pydicom

ROOT = r"E:\NeuroImagepre"

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

def dcm_probe(series_dir):
    """Read the header of the first DICOM file in a series directory"""
    files = sorted(os.listdir(series_dir))
    for fn in files:
        fp = os.path.join(series_dir, fn)
        if not os.path.isfile(fp):
            continue
        try:
            ds = pydicom.dcmread(fp, stop_before_pixels=True, force=True)
        except Exception:
            continue
        if not hasattr(ds, "SOPClassUID") and not hasattr(ds, "Rows"):
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
        # number of files = number of slices
        rec["n_files"] = sum(1 for f in files if os.path.isfile(os.path.join(series_dir, f)))
        return rec
    return None

def main():
    rows = []
    pats = sorted(os.listdir(ROOT))
    for i, pat in enumerate(pats):
        pdir = os.path.join(ROOT, pat)
        if not os.path.isdir(pdir):
            continue
        # study-date directories (there may be several)
        dates = sorted([d for d in os.listdir(pdir) if os.path.isdir(os.path.join(pdir, d))])
        for dt in dates:
            dt_dir = os.path.join(pdir, dt)
            # ImageFileName_* or a direct series directory
            subs = sorted([s for s in os.listdir(dt_dir) if os.path.isdir(os.path.join(dt_dir, s))])
            for sub in subs:
                sub_dir = os.path.join(dt_dir, sub)
                # a numeric directory nested one level deeper means the series level
                inner = sorted([x for x in os.listdir(sub_dir) if os.path.isdir(os.path.join(sub_dir, x))])
                ser_dirs = []
                if inner and any(x.isdigit() for x in inner):
                    ser_dirs = [os.path.join(sub_dir, x) for x in inner if x.isdigit()]
                else:
                    ser_dirs = [sub_dir]
                for sd in ser_dirs:
                    rec = dcm_probe(sd)
                    if rec is None:
                        continue
                    rec["patient"] = pat
                    rec["study_date_dir"] = dt
                    rec["level1"] = sub
                    rec["series_dir"] = os.path.basename(sd)
                    rows.append(rec)
        if (i + 1) % 20 == 0:
            print(f"  ...{i+1}/{len(pats)} patients, {len(rows)} series", flush=True)

    out = r"D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\scripts\pre_dicom_series_inventory.csv"
    cols = ["patient", "study_date_dir", "level1", "series_dir"] + TAGS + ["n_files"]
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nDONE: {len(pats)} patients -> {len(rows)} series -> {out}")

if __name__ == "__main__":
    main()
