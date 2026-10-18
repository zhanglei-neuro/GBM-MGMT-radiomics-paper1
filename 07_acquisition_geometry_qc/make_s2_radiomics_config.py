# -*- coding: utf-8 -*-
r"""
Table S2 — PyRadiomics feature extraction configuration (as actually applied in Paper 1).

配置事实来自两处可核查来源:
  1) 提取脚本  C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts\extract_features.py
  2) 已产出的特征长表 radiomics_features_long.csv（用于反向核对每个类别的特征计数）
"""
import os
import re
import pandas as pd

PROD = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts'
OUT = r'D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1\supplementary\tables'
os.makedirs(OUT, exist_ok=True)

# ---------- 从特征表反向核对各类别计数 ----------
feat = pd.read_csv(os.path.join(PROD, 'radiomics_features_long.csv'))


def cls_of(f):
    if f.startswith('original_shape'):
        return 'shape (3D)'
    m = re.match(r'(original|wavelet-[A-Z]{3})_([a-z0-9]+)_', f)
    return m.group(2) if m else 'other'


def img_of(f):
    if f.startswith('original_shape'):
        return 'shape (3D)'
    m = re.match(r'(original|wavelet-[A-Z]{3})_', f)
    return m.group(1) if m else 'other'


feat['cls'] = feat.feature.map(cls_of)
feat['img'] = feat.feature.map(img_of)
# 注意: nunique 统计的是"全部像型合计"的唯一特征名，需除以波段数还原每像型的类别计数
n_bands = feat[feat.img.str.startswith('wavelet')].img.nunique()          # 8
n_img_types = 1 + n_bands                                                  # 9
cnt_all = feat.groupby('cls').feature.nunique().to_dict()
TEXTURE = ['firstorder', 'glcm', 'glrlm', 'glszm', 'gldm', 'ngtdm']
cnt = {c: cnt_all[c] // n_img_types for c in TEXTURE if c in cnt_all}
cnt['shape'] = cnt_all.get('shape (3D)', 0)
n_basic = sum(cnt[c] for c in TEXTURE)                                     # 每像型 93
per_pair = feat.groupby(['roi', 'seq']).feature.nunique()

assert n_basic == 93, n_basic
assert sum(cnt_all[c] for c in TEXTURE) == n_basic * n_img_types

NAME = {'firstorder': 'First order', 'glcm': 'GLCM (grey level co-occurrence)',
        'glrlm': 'GLRLM (run length)', 'glszm': 'GLSZM (size zone)',
        'gldm': 'GLDM (dependence)', 'ngtdm': 'NGTDM (neighbouring grey tone)',
        'shape': 'Shape (3D)'}
band_order = ['Original'] + ['wavelet-' + b for b in ['LLL', 'LLH', 'LHL', 'LHH', 'HLL', 'HLH', 'HHL', 'HHH']]
# ---------- S2a: 提取设置 ----------
n_img_types = 1 + len([x for x in feat.img.unique() if x.startswith('wavelet')])

s2a_rows = [
    ('Software', 'PyRadiomics 3.1.0', 'Python 3.9 environment (nnUNet virtualenv)'),
    ('Image format', 'NIfTI-1 (.nii.gz)', 'Raidionics pipeline outputs; intensities were NOT modified'),
    ('Mask format', 'NIfTI-1 label map, label value = 1', 'SimpleITK/LabelStatisticsImageFilter'),
    ('Reference grid', 'Contrast-enhanced T1 (T1c)', 'All masks are in the T1c space;'
     ' native T1c is the reference for every patient'),
    ('Non-T1c sequences', 'T1w / FLAIR used on native grid when geometry matched T1c,'
     ' otherwise resampled to T1c grid with ANTs (affine + B-spline warp,'
     ' linear interpolation)', 'ANTs transforms produced by Raidionics'),

    ('Resampling (resampledPixelSpacing)', 'None — no resampling to an isotropic grid',
     'Features are computed at native in-plane resolution and slice thickness'),
    ('Interpolator', 'Not applied (no resampling)', ''),
    ('Intensity normalisation (normalize)', 'None', 'Raw intensities after DICOM→NIfTI conversion'),
    ('Outlier removal (removeOutliers)', 'None', ''),
    ('Voxel array shift (voxelArrayShift)', '0 (default)', ''),
    ('Discretisation', 'Fixed bin width', 'binWidth = 25 (image units)'),
    ('Discretisation (binCount)', 'Not used', 'Fixed bin width takes precedence'),
    ('Resegmentation', 'None', 'Masks used as produced by Raidionics'),
    ('correctMask', 'True', 'Guard for residual sub-voxel geometry differences'),
    ('geometryTolerance', '1e-3', ''),
    ('Voxel-based / force2D', 'False (3D extraction)', 'All ROIs ≥ 30 voxels were retained'),

    ('Enabled image types', 'Original + Wavelet', 'Wavelet = 8 decomposition bands '
     '(LLL, LLH, LHL, LHH, HLL, HLH, HHL, HHH) using PyRadiomics default (coiflet-like) bank'),
    ('Enabled feature classes', 'shape (3D), first order, GLCM, GLRLM, GLSZM, GLDM, NGTDM',
     'shape computed only for the four primary (ROI × sequence) pairings'),
    ('Feature values kept', 'Finite floats only', 'Non-finite values discarded at extraction'),
]
s2a = pd.DataFrame(s2a_rows, columns=['Setting', 'Value', 'Notes'])

# ---------- S2b: 特征类构成 ----------
s2b_rows = [('Shape (3D)', cnt.get('shape', 0), 'original only', 'Primary pairings only')]
for c in ['firstorder', 'glcm', 'glrlm', 'glszm', 'gldm', 'ngtdm']:
    if c in cnt:
        s2b_rows.append((NAME[c], cnt[c], 'original + wavelet', ''))
s2b = pd.DataFrame(s2b_rows, columns=['Feature class', 'N features (per image type)',
                                      'Image types', 'Remarks'])

# ---------- S2c: 像型 × 类别 维度矩阵 ----------
s2c_rows = []
for cls in ['shape', 'firstorder', 'glcm', 'glrlm', 'glszm', 'gldm', 'ngtdm']:
    if cls not in cnt:
        continue
    r = {'Feature class': NAME[cls]}
    r['Original'] = cnt[cls] if cls == 'shape' else cnt[cls]
    for b in band_order[1:]:
        r[b.replace('wavelet-', 'W-')] = '' if cls == 'shape' else cnt[cls]
    r['Row total'] = (cnt[cls] if cls == 'shape' else cnt[cls] * n_img_types)
    s2c_rows.append(r)
s2c = pd.DataFrame(s2c_rows)
tot = pd.DataFrame([{'Feature class': 'TOTAL per (ROI × sequence) pairing',
                     'Original': sum(v for k, v in cnt.items() if k != 'shape'),
                     **{b.replace('wavelet-', 'W-'): sum(v for k, v in cnt.items() if k != 'shape')
                        for b in band_order[1:]},
                     'Row total': sum(v for k, v in cnt.items() if k != 'shape') * n_img_types + cnt.get('shape', 0)}])
s2c = pd.concat([s2c, tot], ignore_index=True)

# ---------- S2d: ROI × 序列配对与总维度 ----------
PRIMARY = [('TumorCore', 'T1c'), ('Necrosis', 'T1c'), ('Edema', 'FLAIR'), ('FLAIRChanges', 'FLAIR')]
CROSS = [('TumorCore', 'FLAIR'), ('TumorCore', 'T1w'), ('Necrosis', 'FLAIR'), ('Necrosis', 'T1w'),
         ('Edema', 'T1c'), ('Edema', 'T1w'), ('FLAIRChanges', 'T1c'), ('FLAIRChanges', 'T1w')]
n_basic = sum(v for k, v in cnt.items() if k != 'shape')          # 93
n_shape = cnt.get('shape', 14)
prim_dim = n_basic * n_img_types + n_shape                        # 851
cross_dim = n_basic                                               # 93
s2d_rows = []
for roi, seq in PRIMARY:
    s2d_rows.append((roi, seq, 'Primary', 'Original + Wavelet (9)', 'shape + 6 texture/intensity classes',
                     prim_dim))
for roi, seq in CROSS:
    s2d_rows.append((roi, seq, 'Cross', 'Original (1)', '6 texture/intensity classes (no shape)',
                     cross_dim))
s2d = pd.DataFrame(s2d_rows, columns=['ROI (Raidionics subregion)', 'Sequence', 'Design',
                                      'Image types', 'Feature classes', 'N features'])
s2d.loc[len(s2d)] = ['— TOTAL —', '4 primary + 8 cross', '', '',
                     'Feature vector per patient', 4 * prim_dim + 8 * cross_dim]

s2a.to_csv(os.path.join(OUT, 'S2a_radiomics_settings.csv'), index=False, encoding='utf-8-sig')
s2b.to_csv(os.path.join(OUT, 'S2b_feature_class_composition.csv'), index=False, encoding='utf-8-sig')
s2c.to_csv(os.path.join(OUT, 'S2c_image_type_matrix.csv'), index=False, encoding='utf-8-sig')
s2d.to_csv(os.path.join(OUT, 'S2d_roi_sequence_pairing.csv'), index=False, encoding='utf-8-sig')

print('=== S2b 特征类构成 ===')
print(s2b.to_string(index=False))
print('\n=== S2d 配对与维度 ===')
print(s2d.to_string(index=False))
print('\n总维度校验: 4 x %d + 8 x %d = %d  (期望 4148)' % (prim_dim, cross_dim, 4 * prim_dim + 8 * cross_dim))
print('每 (roi,seq) 实际计数抽样:')
print(per_pair.head(12).to_string())
