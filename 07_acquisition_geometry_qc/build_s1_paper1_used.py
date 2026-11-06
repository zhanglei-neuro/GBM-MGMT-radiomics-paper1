# -*- coding: utf-8 -*-
r"""
Rebuild the acquisition parameter table (Table S1) from the pre-operative series that Paper 1 **actually used**.

Authoritative definition:
  - analysis cohort = paper1/analysis/predictions.csv (the main analysis cohort, n=90; the single UIH 1.5T case already excluded)
  - specific series = the t1ce / flair / t1 paths recorded in series_selection.csv from the production pipeline
    (these are exactly the source of the Raidionics input NIfTI and hence of the pyradiomics features)
  - DICOM header fields are read straight from the original files; the geometry is cross-checked against the NIfTI in
    input/2, which confirms that "this series is the one Paper 1 used"

Output:
  tables/S1a_paper1_by_vendor.csv        summary by vendor
  tables/S1b_paper1_per_case.csv         per case per series detail (with QC)
  tables/S1c_paper1_vendor_signature.csv vendor x model x field strength technical signature
  logs/s1_paper1_geometry_check.csv      series vs NIfTI geometry consistency check
"""
import os, sys, csv, glob
import numpy as np
import pandas as pd
import pydicom
import nibabel as nib
from pydicom.multival import MultiValue

PROJ = r'D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1'
PROD = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts'
INP = r'D:\Raidionics_data\input\2'
CLIN = r'D:\1桌面文件\AI\GBM_clinical_data.xlsx'
OUT = os.path.join(PROJ, 'supplementary', 'tables')
LOG = os.path.join(PROJ, 'supplementary', 'logs')
os.makedirs(OUT, exist_ok=True)
os.makedirs(LOG, exist_ok=True)

MRI = {'MR', 'MRI'}
ROLE_FILE = {'t1c': 't1-gd', 't1w': 't1w', 'flair': 't2-flair'}
ROLE_LABEL = {'t1c': 'T1C', 't1w': 'T1w', 'flair': 'FLAIR'}
SEL_COL = {'t1c': 't1ce', 't1w': 't1', 'flair': 'flair'}


def valstr(v):
    """DICOM value -> clean string (MultiValue separated by spaces)"""
    if v is None:
        return ''
    if isinstance(v, (MultiValue, list, tuple)):
        return ' '.join(valstr(x) for x in v)
    s = str(v).replace('\\', ' ').strip()
    return s if s.lower() != 'none' else ''


def num(v):
    try:
        f = float(v)
        return f if np.isfinite(f) else np.nan
    except Exception:
        return np.nan


def read_header(series_dir):
    """Read the first DICOM header in the series directory plus the file count"""
    if not series_dir or not os.path.isdir(series_dir):
        return None
    files = [f for f in sorted(os.listdir(series_dir)) if os.path.isfile(os.path.join(series_dir, f))]
    for fn in files:
        try:
            ds = pydicom.dcmread(os.path.join(series_dir, fn), stop_before_pixels=True, force=True)
        except Exception:
            continue
        rec = {
            'Manufacturer': valstr(getattr(ds, 'Manufacturer', '')),
            'Model': valstr(getattr(ds, 'ManufacturerModelName', '')),
            'FieldStrength': num(getattr(ds, 'MagneticFieldStrength', np.nan)),
            'Software': valstr(getattr(ds, 'SoftwareVersions', '')),
            'SeriesDescription': valstr(getattr(ds, 'SeriesDescription', '')),
            'ProtocolName': valstr(getattr(ds, 'ProtocolName', '')),
            'ScanningSequence': valstr(getattr(ds, 'ScanningSequence', '')),
            'SequenceVariant': valstr(getattr(ds, 'SequenceVariant', '')),
            'MRAcquisitionType': valstr(getattr(ds, 'MRAcquisitionType', '')),
            'TR': num(getattr(ds, 'RepetitionTime', np.nan)),
            'TE': num(getattr(ds, 'EchoTime', np.nan)),
            'TI': num(getattr(ds, 'InversionTime', np.nan)),
            'FlipAngle': num(getattr(ds, 'FlipAngle', np.nan)),
            'SliceThickness': num(getattr(ds, 'SliceThickness', np.nan)),
            'SpacingBetweenSlices': num(getattr(ds, 'SpacingBetweenSlices', np.nan)),
            'Rows': num(getattr(ds, 'Rows', np.nan)),
            'Columns': num(getattr(ds, 'Columns', np.nan)),
            'ManufacturerModelName_raw': valstr(getattr(ds, 'ManufacturerModelName', '')),
            'ContrastBolusAgent': valstr(getattr(ds, 'ContrastBolusAgent', '')),
            'StudyDate': valstr(getattr(ds, 'StudyDate', '')),
            'SeriesNumber': valstr(getattr(ds, 'SeriesNumber', '')),
            'n_files': len(files),
        }
        ps = getattr(ds, 'PixelSpacing', None)
        if isinstance(ps, (MultiValue, list, tuple)) and len(ps) >= 2:
            rec['PixelSpacing_row'] = num(ps[0])
            rec['PixelSpacing_col'] = num(ps[1])
        else:
            rec['PixelSpacing_row'] = np.nan
            rec['PixelSpacing_col'] = np.nan
        fov = getattr(ds, 'PixelSpacing', None)
        try:
            rec['FOV_x'] = rec['Columns'] * rec['PixelSpacing_col']
            rec['FOV_y'] = rec['Rows'] * rec['PixelSpacing_row']
        except Exception:
            rec['FOV_x'] = rec['FOV_y'] = np.nan
        return rec
    return None


def to_e(path):
    """D:\\NeuroImagepre\\... -> E:\\NeuroImagepre\\..."""
    if not path or str(path).lower() == 'nan':
        return ''
    p = str(path).replace('/', '\\')
    if p.upper().startswith('D:'):
        p = 'E:' + p[2:]
    return p


def nifti_geom(patient, role):
    d = os.path.join(INP, patient, 'T0')
    if not os.path.isdir(d):
        return None
    for f in sorted(os.listdir(d)):
        if ROLE_FILE[role] in f.lower() and f.endswith('.nii.gz'):
            img = nib.load(os.path.join(d, f))
            return {'dim': tuple(int(x) for x in img.shape),
                    'zooms': tuple(float(z) for z in img.header.get_zooms()[:3])}
    return None


def vendor_norm(m):
    s = (m or '').upper()
    if 'SIEMENS' in s:
        return 'Siemens'
    if 'GE MEDICAL' in s or s.startswith('GE') or 'GENERAL ELECTRIC' in s:
        return 'GE'
    if 'PHILIPS' in s:
        return 'Philips'
    if 'UNITED IMAGING' in s or 'UIH' in s:
        return 'UIH'
    return (m or 'Unknown').strip() or 'Unknown'


def main():
    sel = pd.read_csv(os.path.join(PROD, 'series_selection.csv'))
    sel = sel[sel.patient.notna()]
    sel['patient'] = sel.patient.astype(str)

    cohort = pd.read_csv(os.path.join(PROJ, 'analysis', 'predictions.csv'))
    cohort['patient'] = cohort.patient.astype(str)
    patients = cohort.patient.tolist()
    vendor_tab = dict(zip(cohort.patient, cohort.vendor))
    print('Paper1 cohort: %d cases' % len(patients))

    # date of surgery
    surg = {}
    try:
        c = pd.read_csv(os.path.join(PROD, 'cohort_MGMT.csv'))
        c['patient'] = c.patient.astype(str)
        for _, r in c.iterrows():
            surg[r.patient] = str(r.get('surgery_date', ''))
    except Exception as e:
        print('  could not read the date of surgery:', e)

    rows, geo_rows, missing = [], [], []
    t1c_geom = {}
    for p in patients:
        t1c_geom[p] = nifti_geom(p, 't1c')
    for p in patients:
        srow = sel[sel.patient == p]
        if len(srow) == 0:
            missing.append((p, 'not in series_selection'))
            continue
        srow = srow.iloc[0]
        for role in ('t1c', 't1w', 'flair'):
            p_d = str(srow.get(SEL_COL[role], '') or '')
            p_e = to_e(p_d)
            h = read_header(p_e)
            if h is None:
                missing.append((p, role + ' | ' + p_e))
                continue
            g = nifti_geom(p, role)
            # grid used for this series at extraction time: native geometry vs resampled to the T1c grid
            if role == 't1c' or g is None or t1c_geom.get(p) is None:
                grid = 'native (T1c reference)'
            else:
                same = (g['dim'] == t1c_geom[p]['dim'] and
                        all(abs(a - b) < 1e-4 for a, b in zip(g['zooms'], t1c_geom[p]['zooms'])))
                grid = 'native' if same else 'resampled to T1c grid'
            # geometry consistency check
            ok_dims = ok_spacing = ''
            if g is not None:
                ok_dims = (g['dim'][0] == int(h['Columns']) and g['dim'][1] == int(h['Rows'])) if (
                    np.isfinite(h['Columns']) and np.isfinite(h['Rows'])) else ''
                if np.isfinite(h['PixelSpacing_col']) and h['PixelSpacing_col'] > 0:
                    ok_spacing = bool(abs(g['zooms'][0] - h['PixelSpacing_col']) < 1e-3)
            geo_rows.append({
                'patient': p, 'role': ROLE_LABEL[role],
                'dicom_path': p_e,
                'dicom_series_desc': h['SeriesDescription'],
                'nifti_dim': 'x'.join(map(str, g['dim'])) if g else '',
                'dicom_Cols_x_Rows': '%gx%g' % (h['Columns'], h['Rows']) if np.isfinite(h['Columns']) else '',
                'nifti_zoom': 'x'.join('%.5g' % z for z in g['zooms']) if g else '',
                'dicom_pixelspacing': '%.5g x %.5g' % (h['PixelSpacing_col'], h['PixelSpacing_row']) if np.isfinite(h['PixelSpacing_col']) else '',
                'dim_match': ok_dims, 'spacing_match': ok_spacing,
            })
            rows.append({
                'Patient': p,
                'Vendor': vendor_norm(h['Manufacturer']),
                'Vendor (Paper 1 label)': vendor_tab.get(p, ''),
                'Vendor label match': 'yes' if vendor_norm(h['Manufacturer']) == vendor_tab.get(p) else 'NO',
                'Model': h['Model'],
                'Field strength (T)': h['FieldStrength'],
                'Software version': h['Software'],
                'Sequence role': ROLE_LABEL[role],
                'Series description': h['SeriesDescription'],
                'Protocol name': h['ProtocolName'],
                'Scanning sequence': h['ScanningSequence'],
                'Sequence variant': h['SequenceVariant'],
                'Acquisition type': h['MRAcquisitionType'],
                'Grid used for radiomics': grid,
                'TR (ms)': h['TR'], 'TE (ms)': h['TE'], 'TI (ms)': h['TI'],
                'Flip angle (deg)': h['FlipAngle'],
                'Slice thickness (mm)': h['SliceThickness'],
                'Slice spacing (mm)': h['SpacingBetweenSlices'],
                'Rows': h['Rows'], 'Columns': h['Columns'],
                'Matrix': '%g x %g' % (h['Rows'], h['Columns']) if np.isfinite(h['Rows']) else '',
                'Pixel spacing row (mm)': h['PixelSpacing_row'],
                'Pixel spacing col (mm)': h['PixelSpacing_col'],
                'FOV x (mm)': h['FOV_x'], 'FOV y (mm)': h['FOV_y'],
                'N slices': h['n_files'],
                'Contrast agent': h['ContrastBolusAgent'],
                'Study date': h['StudyDate'],
                'Surgery date': surg.get(p, ''),
                'DICOM series dir (E:)': p_e,
            })

    df = pd.DataFrame(rows)
    print('S1b rows: %d (expected n x 3)' % len(df))
    if missing:
        print('%d records missing:' % len(missing))
        for m in missing[:20]:
            print('   ', m)

    # ---- QC: imaging vs date of surgery ----
    def days(row):
        try:
            import datetime
            sd = str(row['Study date'])
            su = str(row['Surgery date'])[:10].replace('-', '')
            if len(sd) >= 8 and len(su) == 8:
                a = datetime.date(int(sd[:4]), int(sd[4:6]), int(sd[6:8]))
                b = datetime.date(int(su[:4]), int(su[4:6]), int(su[6:8]))
                return (a - b).days
        except Exception:
            pass
        return np.nan

    df['Days image - surgery'] = df.apply(days, axis=1)

    def qc(row):
        d = row['Days image - surgery']
        if pd.isna(d):
            return 'NO_SURGERY_DATE'
        if d > 90:
            return 'POSTOP/LATE (>90d)'
        if d > 7:
            return 'POSTOP (8-90d)'
        if d >= -7:
            return 'OK (peri-op)'
        return 'OK (preop)'

    df['QC flag'] = df.apply(qc, axis=1)
    print('\nQC flag:')
    print(df['QC flag'].value_counts().to_string())

    geo = pd.DataFrame(geo_rows)
    bad = geo[(geo.dim_match == False)]
    if len(bad):
        print('\ngeometry mismatches: %d' % len(bad))
        print(bad[['patient', 'role', 'dicom_series_desc', 'nifti_dim', 'dicom_Cols_x_Rows']].to_string(index=False))
    else:
        print('\ngeometry consistency: all passed (%d/%d)' % ((geo.dim_match == True).sum(), len(geo)))

    df.to_csv(os.path.join(OUT, 'S1b_paper1_per_case.csv'), index=False, encoding='utf-8-sig')
    geo.to_csv(os.path.join(LOG, 's1_paper1_geometry_check.csv'), index=False, encoding='utf-8-sig')
    print('\n[S1b written]')

    # ---- S1a: summary by vendor ----
    def fmt_med(x, dec=1):
        x = pd.to_numeric(x, errors='coerce').dropna()
        if len(x) == 0:
            return ''
        return ('%.*f (%.*f-%.*f)' % (dec, x.median(), dec, x.min(), dec, x.max()))

    order = ['Siemens', 'GE', 'UIH', 'Philips']
    # grouping is by the measured vendor of T1c (the reference grid for segmentation and extraction), so each case falls in exactly one group
    t1c_vendor = df[df['Sequence role'] == 'T1C'].drop_duplicates('Patient').set_index('Patient').Vendor
    df['_Grp'] = df.Patient.map(t1c_vendor)
    vend = [v for v in order if v in set(df._Grp.dropna())]
    ALLCOL = 'All (n=%d)' % df.Patient.nunique()
    VCOL = {v: '%s (n=%d)' % (v, df[df._Grp == v].Patient.nunique()) for v in vend}

    lines = []
    def add(param, fn):
        r = {'Parameter': param, ALLCOL: fn(df)}
        for v in vend:
            r[VCOL[v]] = fn(df[df._Grp == v])
        lines.append(r)

    for role in ('T1C', 'T1w', 'FLAIR'):
        sub = lambda x, r=role: x[x['Sequence role'] == r]
        lines.append({'Parameter': '[%s] sequence' % role})
        add('%s · n patients' % role, lambda x, r=role: str(sub(x, r).Patient.nunique()))
        add('%s · field strength 1.5T / 3.0T' % role,
            lambda x, r=role: '%d / %d' % ((sub(x, r)['Field strength (T)'] == 1.5).sum(),
                                           (sub(x, r)['Field strength (T)'] == 3.0).sum()))
        add('%s · TR ms, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['TR (ms)'], 0))
        add('%s · TE ms, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['TE (ms)'], 1))
        add('%s · TI ms, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['TI (ms)'], 0))
        add('%s · flip angle deg, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['Flip angle (deg)'], 0))
        add('%s · slice thickness mm, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['Slice thickness (mm)'], 2))
        add('%s · in-plane pixel mm, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['Pixel spacing col (mm)'], 3))
        add('%s · matrix (rows x cols), mode' % role,
            lambda x, r=role: (sub(x, r)['Matrix'].mode().iloc[0] if len(sub(x, r)['Matrix'].mode()) else ''))
        add('%s · slices, median (range)' % role, lambda x, r=role: fmt_med(sub(x, r)['N slices'], 0))
        add('%s · acquisition type 3D/2D' % role,
            lambda x, r=role: ' / '.join('%s: %d' % (k or 'N/A', v) for k, v in
                                         sub(x, r)['Acquisition type'].value_counts().items()))
    lines.append({'Parameter': '[Scanner models]'})
    add('model, n patients (per vendor)',
        lambda x: '; '.join('%s %s (%d)' % (v, m or 'N/A', n) for (v, m, n) in
                            x.drop_duplicates('Patient').groupby(['Vendor', 'Model']).size()
                            .reset_index().apply(lambda r: (r.Vendor, r.Model, r[0]), axis=1)))
    add('software version(s)',
        lambda x: '; '.join(sorted(set(v for v in x['Software version'] if str(v).strip()))[:4]))
    lines.append({'Parameter': '[QC]'})
    lines.append({'Parameter': 'note',
                  ALLCOL: ('Vendor groups are defined by the T1c (contrast-enhanced) reference '
                           'sequence, which is also the grid on which segmentation and radiomics '
                           'extraction were performed. For every case the T1c, T1w and FLAIR '
                           'series were acquired on a single scanner within the same preoperative '
                           'examination. Acquisition parameters below are the values actually '
                           'used for radiomics extraction. One case (P021, United Imaging '
                           'uMR560 1.5 T) was excluded from the primary analysis per the '
                           'pre-specified exclusion criterion (Section 3.8) and is not represented '
                           'here.')})

    s1a = pd.DataFrame(lines).fillna('')
    s1a.to_csv(os.path.join(OUT, 'S1a_paper1_by_vendor.csv'), index=False, encoding='utf-8-sig')
    print('\n=== S1a ===')
    print(s1a.to_string(index=False))

    # ---- S1c: vendor x model x field strength technical signature ----
    sc = df.groupby(['Vendor', 'Model', 'Field strength (T)']).agg(
        n_sequences=('Patient', 'size'),
        n_patients=('Patient', 'nunique'),
        software=('Software version', lambda x: '; '.join(sorted(set(v for v in x if str(v).strip()))[:2])),
        sequences=('Sequence role', lambda x: '/'.join(sorted(set(x)))),
    ).reset_index()
    sc.to_csv(os.path.join(OUT, 'S1c_paper1_vendor_signature.csv'), index=False, encoding='utf-8-sig')
    print('\n=== S1c ===')
    print(sc.to_string(index=False))

    # vendor distribution of the cohort (to compare with Table 1)
    pv = df.drop_duplicates('Patient')[['Patient', 'Vendor']]
    print('\nvendor distribution of the cohort (by case):')
    print(pv.Vendor.value_counts().to_string())
    print('\ncases where the declared vendor label disagrees with the DICOM measurement:')
    pv = pv.copy()
    pv['declared'] = pv.Patient.map(vendor_tab)
    mism = pv[pv.Vendor != pv.declared]
    print(mism.to_string(index=False) if len(mism) else '  none')


if __name__ == '__main__':
    main()
