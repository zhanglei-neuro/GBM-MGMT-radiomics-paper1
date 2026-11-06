# -*- coding: utf-8 -*-
r"""Raidionics GBM batch segmentation - step 1: survey the DICOM data

Scan every patient series under D:\NeuroImagepre and produce:
  - dicom_inventory.csv   header information for each series
  - series_selection.csv  the T1-CE / FLAIR / unenhanced T1 series path chosen for each patient

Decision rules (hard rules):
  - T1-CE is decided by a non-empty ContrastBolusAgent (0018,0010); never trust the series name
  - missing-slice detection uses the projected spacing along the slice normal, not InstanceNumber or a
    simple z difference
"""
import os, sys, csv, math, traceback
import pydicom
import numpy as np

ROOT = r"D:\NeuroImagepre"
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

def find_series_dirs(patient_dir):
    """Return [(series_dir, [dcm_files])]"""
    result = []
    for dirpath, dirnames, filenames in os.walk(patient_dir):
        dcms = [f for f in filenames if f.lower().endswith(('.dcm', '.ima'))]
        if dcms:
            result.append((dirpath, sorted(dcms)))
    return result

def series_info(series_dir, files):
    """Read the series header information"""
    n = len(files)
    first = pydicom.dcmread(os.path.join(series_dir, files[0]), stop_before_pixels=True)
    mid = pydicom.dcmread(os.path.join(series_dir, files[n // 2]), stop_before_pixels=True)

    desc = str(getattr(first, 'SeriesDescription', '') or '')
    series_num = str(getattr(first, 'SeriesNumber', '') or '')
    study_date = str(getattr(first, 'StudyDate', '') or '')
    modality = str(getattr(first, 'Modality', '') or '')
    rows = int(getattr(first, 'Rows', 0) or 0)
    cols = int(getattr(first, 'Columns', 0) or 0)
    agent = str(getattr(first, 'ContrastBolusAgent', '') or '').strip()

    iop = getattr(first, 'ImageOrientationPatient', None)
    ipp_first = getattr(first, 'ImagePositionPatient', None)
    ipp_mid = getattr(mid, 'ImagePositionPatient', None)

    px_sp = None
    try:
        ps = first.PixelSpacing
        px_sp = (float(ps[0]), float(ps[1]))
    except Exception:
        pass

    slice_sp = None
    gap_bad = False
    gap_report = ''
    try:
        if iop is not None and ipp_first is not None and ipp_mid is not None:
            normal = np.cross(np.array(iop[:3], dtype=float), np.array(iop[3:], dtype=float))
            normal /= np.linalg.norm(normal)
            # sample projections over several slices
            positions = []
            sample_idx = list(range(0, n, max(1, n // 40)))
            for i in sample_idx:
                d = pydicom.dcmread(os.path.join(series_dir, files[i]), stop_before_pixels=True)
                p = getattr(d, 'ImagePositionPatient', None)
                if p is None:
                    continue
                positions.append(float(np.dot(np.array([float(p[0]), float(p[1]), float(p[2])]), normal)))
            positions.sort()
            if len(positions) >= 3:
                diffs = np.diff(positions)
                diffs = diffs[diffs > 0.01]  # collapse repeated positions
                if len(diffs):
                    slice_sp = round(float(np.median(diffs)), 3)
                    # missing-slice detection: spacing variation > 30%
                    if diffs.max() > slice_sp * 1.35 or diffs.min() < slice_sp * 0.65:
                        gap_bad = True
                        gap_report = f"min={diffs.min():.2f} max={diffs.max():.2f} med={slice_sp}"
        # without an orientation tag, fall back to the z difference
        if slice_sp is None and ipp_first is not None and ipp_mid is not None:
            dz = abs(float(ipp_mid[2]) - float(ipp_first[2]))
            slice_sp = round(dz, 3) if dz > 0.01 else None
    except Exception as e:
        gap_report = 'ERR:' + str(e)[:60]

    return dict(series_dir=series_dir, n_files=n, series_num=series_num, desc=desc,
                study_date=study_date, modality=modality, rows=rows, cols=cols,
                contrast_agent=agent, px_sp=px_sp, slice_sp=slice_sp,
                gap_bad=gap_bad, gap_report=gap_report)

def classify(info):
    """Role assignment: t1ce / flair / t1 / other"""
    desc_l = info['desc'].lower()
    has_contrast = bool(info['contrast_agent'])
    is_t1 = ('t1' in desc_l) or ('mprage' in desc_l) or ('mp-rage' in desc_l)
    is_flair = ('flair' in desc_l) and ('t1' not in desc_l)  # "T1 FLAIR" is often a decoy name for an unenhanced series
    if not is_flair and ('flair' in desc_l):
        # where the name holds both t1 and flair, let the contrast tag decide
        is_flair = False
    if is_t1 and has_contrast:
        return 't1ce'
    if is_flair and not has_contrast:
        return 'flair'
    if is_t1 and not has_contrast:
        return 't1'
    if has_contrast and not is_flair:
        return 't1ce?'  # contrast present but no T1 in the name: probably T1-CE, needs a manual check
    return 'other'

def main():
    patients = sorted([d for d in os.listdir(ROOT)
                       if os.path.isdir(os.path.join(ROOT, d))])
    inv_rows, sel_rows = [], []
    for pi, pat in enumerate(patients):
        try:
            series_list = find_series_dirs(os.path.join(ROOT, pat))
        except Exception as e:
            print(f"[ERR] {pat}: {e}")
            continue
        infos = []
        for sdir, files in series_list:
            try:
                info = series_info(sdir, files)
            except Exception as e:
                print(f"[ERR] {pat} {sdir}: {e}")
                continue
            info['patient'] = pat
            info['role'] = classify(info)
            infos.append(info)
            inv_rows.append(info)
        # selection: per role take the one with most slices and no missing slices; for t1ce prefer the contrast-enhanced one
        def pick(role):
            cand = [x for x in infos if x['role'] == role or (role == 't1ce' and x['role'] == 't1ce?')]
            if not cand:
                return None
            cand.sort(key=lambda x: (x['gap_bad'], -x['n_files']))
            return cand[0]
        sel = {'patient': pat}
        for role in ('t1ce', 'flair', 't1'):
            c = pick(role)
            sel[role] = c['series_dir'] if c else ''
            sel[role + '_desc'] = (c['desc'] + f" ({c['n_files']}slices)") if c else ''
            sel[role + '_gap'] = ('BAD:' + c['gap_report']) if (c and c['gap_bad']) else ('ok' if c else '-')
        sel_rows.append(sel)
        print(f"[{pi+1}/{len(patients)}] {pat}: t1ce={'Y' if sel['t1ce'] else 'N'} "
              f"flair={'Y' if sel['flair'] else 'N'} t1={'Y' if sel['t1'] else 'N'}")

    with open(os.path.join(OUT_DIR, 'dicom_inventory.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['patient', 'series_dir', 'series_num', 'desc', 'study_date',
                                          'modality', 'n_files', 'rows', 'cols', 'contrast_agent',
                                          'px_sp', 'slice_sp', 'role', 'gap_bad', 'gap_report'])
        w.writeheader()
        for r in inv_rows:
            r2 = dict(r)
            r2['px_sp'] = str(r2['px_sp'] or '')
            w.writerow(r2)

    with open(os.path.join(OUT_DIR, 'series_selection.csv'), 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['patient', 't1ce', 't1ce_desc', 't1ce_gap',
                                          'flair', 'flair_desc', 'flair_gap',
                                          't1', 't1_desc', 't1_gap'])
        w.writeheader()
        w.writerows(sel_rows)

    n_complete = sum(1 for s in sel_rows if s['t1ce'] and s['flair'] and s['t1'])
    print(f"\n{len(patients)} patients in total, {n_complete} with all three sequences")
    print("done")

if __name__ == '__main__':
    main()
