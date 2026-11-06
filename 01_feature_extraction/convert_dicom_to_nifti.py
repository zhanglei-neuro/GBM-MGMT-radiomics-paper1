# -*- coding: utf-8 -*-
r"""Steps 2+3: rebuild the input directory and convert DICOM -> NIfTI (file names carry the role keyword)

Input:  series_selection.csv
Output: D:\Raidionics_data\input\2\<Patient>\T0\
  - <pat>_t1-gd_pre.nii.gz    (T1-CE; 'gd' in the name tells the backend this is T1c)
  - <pat>_t2-flair_pre.nii.gz (FLAIR; 'flair' in the name tells the backend this is FLAIR)
  - <pat>_t1w_pre.nii.gz      (unenhanced T1; no keyword, so T1w)

Copy only: read the original DICOM, never resample, normalise or skull-strip.
Resumable: existing target files are skipped.
"""
import csv, os, sys
import SimpleITK as sitk

BASE = os.path.dirname(os.path.abspath(__file__))
SEL = os.path.join(BASE, 'series_selection.csv')
CONV_LOG = os.path.join(BASE, 'conversion_log.csv')
OUT_ROOT = r"D:\Raidionics_data\input\2"

FORBIDDEN = ['gd', 'flair', 't2', 'dwi', 'tirm']  # patient names must not contain these words, or the file names become ambiguous

def resolve_dicom(path):
    """Drive fallback: the original DICOM has moved from D: to E:, but series_selection.csv keeps the historic D: paths.
    Use the original path when it is readable, otherwise try the same path on E:."""
    if not path:
        return path
    if os.path.isdir(path):
        return path
    p = str(path).replace('/', '\\')
    if p.upper().startswith('D:'):
        alt = 'E:' + p[2:]
        if os.path.isdir(alt):
            return alt
    return path


def manual_assemble(series_dir, files_sorted, out_path):
    """Fallback: when slice-position jitter exceeds the ITK tolerance, assemble the volume manually from the DICOM tags (a standard reconstruction, not a resampling)"""
    import pydicom
    import numpy as np
    ds0 = pydicom.dcmread(os.path.join(series_dir, files_sorted[0]))
    iop = np.array(ds0.ImageOrientationPatient, dtype=float)
    normal = np.cross(iop[:3], iop[3:])
    normal /= np.linalg.norm(normal)
    rows, cols = int(ds0.Rows), int(ds0.Columns)
    ps_r, ps_c = float(ds0.PixelSpacing[0]), float(ds0.PixelSpacing[1])
    slices, projections = [], []
    for f in files_sorted:
        ds = pydicom.dcmread(os.path.join(series_dir, f))
        arr = ds.pixel_array.astype(np.float32)
        arr = arr * float(getattr(ds, 'RescaleSlope', 1)) + float(getattr(ds, 'RescaleIntercept', 0))
        slices.append(arr)
        p = np.array(ds.ImagePositionPatient, dtype=float)
        projections.append(float(np.dot(p, normal)))
    proj = np.array(projections)
    diffs = np.diff(proj)
    z_sp = float(np.median(diffs))
    spacing = (ps_c, ps_r, z_sp)
    origin = tuple(float(v) for v in np.array(ds0.ImagePositionPatient, dtype=float))
    # direction matrix: columns are the row direction, the column direction and the normal (LPS)
    direction = tuple(np.stack([iop[3:], iop[:3], normal], axis=1).flatten().tolist())
    vol = np.stack(slices, axis=0)
    img = sitk.GetImageFromArray(vol)
    img.SetSpacing(spacing)
    img.SetOrigin(origin)
    img.SetDirection(direction)
    sitk.WriteImage(img, out_path, True)
    return img.GetSize(), img.GetSpacing()

def convert_series(series_dir, out_path):
    import pydicom
    import numpy as np
    # sort by DICOM geometry and de-duplicate by slice position, so a repeated-slice series cannot pollute the volume
    files = [f for f in os.listdir(series_dir) if f.lower().endswith(('.dcm', '.ima'))]
    if not files:
        raise ValueError('no DICOM files found')
    first = pydicom.dcmread(os.path.join(series_dir, files[0]), stop_before_pixels=True)
    iop = np.array(first.ImageOrientationPatient, dtype=float)
    normal = np.cross(iop[:3], iop[3:])
    normal /= np.linalg.norm(normal)
    entries = []
    for f in files:
        d = pydicom.dcmread(os.path.join(series_dir, f), stop_before_pixels=True)
        p = getattr(d, 'ImagePositionPatient', None)
        if p is None:
            raise ValueError(f'missing ImagePositionPatient in {f}')
        key = (round(float(p[0]), 2), round(float(p[1]), 2), round(float(p[2]), 2))
        shape = (int(d.Rows), int(d.Columns))
        entries.append((float(np.dot(np.array([float(p[0]), float(p[1]), float(p[2])]), normal)), key, f, shape))
    # drop odd slices whose matrix size differs from the mode (screenshots and the like added by the export tool)
    from collections import Counter
    shape_counts = Counter(e[3] for e in entries)
    major_shape = shape_counts.most_common(1)[0][0]
    entries = [e for e in entries if e[3] == major_shape]
    entries.sort(key=lambda x: x[0])
    seen = set()
    unique_files = []
    for _, key, f, _shape in entries:
        if key not in seen:
            seen.add(key)
            unique_files.append(f)
    fl = [os.path.join(series_dir, f) for f in unique_files]
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(fl)
    try:
        img = reader.Execute()
        mode = 'sitk'
    except Exception:
        # ITK tolerates no slice non-uniformity -> assemble manually
        img = None
        mode = 'manual'
    if img is None:
        size, sp = manual_assemble(series_dir, unique_files, out_path)
        return size, sp, len(files), len(unique_files), mode
    sitk.WriteImage(img, out_path, True)
    return img.GetSize(), img.GetSpacing(), len(files), len(unique_files), mode

def main():
    sel = list(csv.DictReader(open(SEL, encoding='utf-8-sig')))
    log_rows = []
    for s in sel:
        p = s['patient']
        pl = p.lower()
        if any(t in pl for t in FORBIDDEN):
            log_rows.append({'patient': p, 'status': 'SKIP_NAME_CONFLICT', 'note': f"patient name contains keyword {[t for t in FORBIDDEN if t in pl]}"})
            print(f"[SKIP] {p}: name conflict")
            continue
        if not (s['t1ce'] and s['flair'] and s['t1']):
            log_rows.append({'patient': p, 'status': 'INCOMPLETE', 'note': 'missing series'})
            print(f"[SKIP] {p}: incomplete")
            continue
        t0 = os.path.join(OUT_ROOT, p, 'T0')
        os.makedirs(t0, exist_ok=True)
        targets = [
            ('t1ce', s['t1ce'], os.path.join(t0, f'{pl}_t1-gd_pre.nii.gz')),
            ('flair', s['flair'], os.path.join(t0, f'{pl}_t2-flair_pre.nii.gz')),
            ('t1', s['t1'], os.path.join(t0, f'{pl}_t1w_pre.nii.gz')),
        ]
        ok = True
        note_parts = []
        for role, src, dst in targets:
            if os.path.exists(dst) and os.path.getsize(dst) > 1000:
                note_parts.append(f'{role}=cached')
                continue
            try:
                src = resolve_dicom(src)
                size, sp, n_raw, n_uni, mode = convert_series(src, dst)
                dup = f' DEDUP({n_raw}->{n_uni})' if n_raw != n_uni else ''
                md = ' MANUAL' if mode == 'manual' else ''
                note_parts.append(f'{role}={size}@{round(sp[0],2)}x{round(sp[2],2)}{dup}{md}')
            except Exception as e:
                ok = False
                log_rows.append({'patient': p, 'status': 'ERROR', 'note': f'{role}: {str(e)[:120]}'})
                print(f"[ERR] {p} {role}: {str(e)[:120]}")
                break
        if ok:
            # cross-check the three volumes for consistent size, spacing and orientation
            try:
                imgs = [sitk.ReadImage(dst) for _, _, dst in targets]
                sizes = {i.GetSize() for i in imgs}
                spaces = {tuple(round(v, 2) for v in i.GetSpacing()) for i in imgs}
                dirs = {tuple(round(v, 4) for v in i.GetDirection()) for i in imgs}
                n_diff = len(sizes) + len(spaces) + len(dirs) - 3
                log_rows.append({'patient': p, 'status': 'OK', 'note': '; '.join(note_parts) +
                                 (f'; MISMATCH={n_diff}' if n_diff else '')})
                print(f"[OK] {p}: {'; '.join(note_parts)}")
            except Exception as e:
                log_rows.append({'patient': p, 'status': 'ERROR', 'note': f'verify: {str(e)[:120]}'})
                print(f"[ERR] {p} verify: {str(e)[:120]}")
        with open(CONV_LOG, 'w', newline='', encoding='utf-8-sig') as f:
            w = csv.DictWriter(f, fieldnames=['patient', 'status', 'note'])
            w.writeheader()
            w.writerows(log_rows)

if __name__ == '__main__':
    main()
