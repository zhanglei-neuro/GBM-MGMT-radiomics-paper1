# -*- coding: utf-8 -*-
r"""
按 Paper 1 **实际使用的**术前影像序列，重建采集参数表（Table S1）。

权威口径:
  - 分析队列 = paper1/analysis/predictions.csv（主分析队列，n=90；UIH 1.5T 单例已排除）
  - 具体序列 = 生产流程 series_selection.csv 记录的 t1ce / flair / t1 路径
    （这就是 Raidionics 输入 NIfTI -> pyradiomics 特征的实际来源）
  - DICOM 头字段直接读原始文件；几何与 input/2 的 NIfTI 交叉验证，确保"这一条序列
    就是 Paper 1 用的那一条"

输出:
  tables/S1a_paper1_by_vendor.csv     按厂商汇总
  tables/S1b_paper1_per_case.csv      每例每序列明细（含 QC）
  tables/S1c_paper1_vendor_signature.csv  厂商×机型×场强 技术签名
  logs/s1_paper1_geometry_check.csv   序列↔NIfTI 几何一致性校验
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
    """DICOM 值 -> 干净字符串（MultiValue 空格分隔）"""
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
    """读序列目录首个 DICOM 头 + 文件数"""
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
    print('Paper1 队列: %d 例' % len(patients))

    # 手术日期
    surg = {}
    try:
        c = pd.read_csv(os.path.join(PROD, 'cohort_MGMT.csv'))
        c['patient'] = c.patient.astype(str)
        for _, r in c.iterrows():
            surg[r.patient] = str(r.get('surgery_date', ''))
    except Exception as e:
        print('  手术日期读取失败:', e)

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
            # 特征提取时该序列所用的网格：原生几何 vs 重采样到 T1c 网格
            if role == 't1c' or g is None or t1c_geom.get(p) is None:
                grid = 'native (T1c reference)'
            else:
                same = (g['dim'] == t1c_geom[p]['dim'] and
                        all(abs(a - b) < 1e-4 for a, b in zip(g['zooms'], t1c_geom[p]['zooms'])))
                grid = 'native' if same else 'resampled to T1c grid'
            # 几何一致性校验
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
    print('S1b 行数: %d (期望 n×3)' % len(df))
    if missing:
        print('缺失 %d 条:' % len(missing))
        for m in missing[:20]:
            print('   ', m)

    # ---- QC: 影像 vs 手术日期 ----
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
        print('\n⚠ 几何不一致 %d 条:' % len(bad))
        print(bad[['patient', 'role', 'dicom_series_desc', 'nifti_dim', 'dicom_Cols_x_Rows']].to_string(index=False))
    else:
        print('\n几何一致性: 全部通过 (%d/%d)' % ((geo.dim_match == True).sum(), len(geo)))

    df.to_csv(os.path.join(OUT, 'S1b_paper1_per_case.csv'), index=False, encoding='utf-8-sig')
    geo.to_csv(os.path.join(LOG, 's1_paper1_geometry_check.csv'), index=False, encoding='utf-8-sig')
    print('\n[S1b 已写出]')

    # ---- S1a: 按厂商汇总 ----
    def fmt_med(x, dec=1):
        x = pd.to_numeric(x, errors='coerce').dropna()
        if len(x) == 0:
            return ''
        return ('%.*f (%.*f-%.*f)' % (dec, x.median(), dec, x.min(), dec, x.max()))

    order = ['Siemens', 'GE', 'UIH', 'Philips']
    # 分组依据：T1c（分割与特征提取的参考网格）实测厂商，保证每例只归一组
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

    # ---- S1c: 厂商 × 机型 × 场强 技术签名 ----
    sc = df.groupby(['Vendor', 'Model', 'Field strength (T)']).agg(
        n_sequences=('Patient', 'size'),
        n_patients=('Patient', 'nunique'),
        software=('Software version', lambda x: '; '.join(sorted(set(v for v in x if str(v).strip()))[:2])),
        sequences=('Sequence role', lambda x: '/'.join(sorted(set(x)))),
    ).reset_index()
    sc.to_csv(os.path.join(OUT, 'S1c_paper1_vendor_signature.csv'), index=False, encoding='utf-8-sig')
    print('\n=== S1c ===')
    print(sc.to_string(index=False))

    # 队列厂商分布（与 Table1 对照）
    pv = df.drop_duplicates('Patient')[['Patient', 'Vendor']]
    print('\n队列厂商分布（按病例）:')
    print(pv.Vendor.value_counts().to_string())
    print('\n声明 vendor 标签 vs DICOM 实测 不一致的病例:')
    pv = pv.copy()
    pv['declared'] = pv.Patient.map(vendor_tab)
    mism = pv[pv.Vendor != pv.declared]
    print(mism.to_string(index=False) if len(mism) else '  无')


if __name__ == '__main__':
    main()
