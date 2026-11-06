# -*- coding: utf-8 -*-
r"""Vendor and field-strength audit for the whole cohort: read the first readable DICOM header under
E:\NeuroImagepre\<patient>\ for every case, and report Manufacturer / ManufacturerModelName / MagneticFieldStrength.

Background: the MRI_vendor column of GBM_clinical_data.xlsx carries stale or wrong labels (P055 is labelled UIH but is
in fact Siemens), so it cannot be relied on for platform attribution, and Paper 1 now uses the DICOM-measured vendor.
This script generalises that practice to every patient and turned up the second United Imaging uMR560 1.5 T case
(P055) hiding in the Paper 2 follow-up cohort.

Usage: python audit_vendor_fieldstrength.py <patient list csv> <output csv>
The patient list csv must contain a column named 'patient'.
"""
import os
import sys
import pandas as pd
import pydicom

DICOM_ROOT = r'E:\NeuroImagepre'


def scan_first_hdr(base):
    """Walk the directory and return (Manufacturer, Model, FieldStrength) from the first readable DICOM header."""
    for root, _dirs, files in os.walk(base):
        for f in sorted(files):
            try:
                ds = pydicom.dcmread(os.path.join(root, f), stop_before_pixels=True)
                man = str(getattr(ds, 'Manufacturer', ''))
                if man:
                    return man, str(getattr(ds, 'ManufacturerModelName', '')), \
                        str(getattr(ds, 'MagneticFieldStrength', ''))
            except Exception:
                pass
    return '', '', ''


def main():
    ids_csv, out_csv = sys.argv[1], sys.argv[2]
    pids = sorted(set(pd.read_csv(ids_csv)['patient'].astype(str)))
    rows = []
    for p in pids:
        base = os.path.join(DICOM_ROOT, p)
        if not os.path.isdir(base):
            rows.append((p, 'NO_DIR', '', ''))
            continue
        rows.append((p,) + scan_first_hdr(base))
    df = pd.DataFrame(rows, columns=['patient', 'Manufacturer', 'Model', 'FieldStrength'])
    df.to_csv(out_csv, index=False, encoding='utf-8-sig')
    print(df.groupby(['Manufacturer', 'Model', 'FieldStrength']).size().to_string())
    bad = df[~df['FieldStrength'].astype(str).str.startswith('3')]
    if len(bad):
        print('\ncases that are not 3.0 T (to be excluded):')
        print(bad.to_string())


if __name__ == '__main__':
    main()
