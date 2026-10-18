# -*- coding: utf-8 -*-
r"""阶段1: pyradiomics 特征提取 (在 D:\Anaconda\envs\nnunet 环境运行, python 3.9 + pyradiomics 3.1.0)

设计:
- 全部掩码都在 t1-gd 空间 (Raidionics 输出, 已验证)
- 序列对齐: t1-gd 用原始; t1w/t2-flair 若与 t1-gd 几何完全一致则用原始,
  否则用 Raidionics 自带的 ANTs 前向变换重采样到 t1-gd 网格
  (CompositeTransform([affine, warp]), 已用 MI 校验为最优组合), 结果缓存到 input_reg
- 掩码用 nibabel 读数组(绕开 SimpleITK 2.5.2 读标注文件崩溃的 bug), 几何从参考影像拷贝
- 配对:
  primary(Original+Wavelet+shape): TumorCore×T1CE, Necrosis×T1CE, Edema×FLAIR, FLAIRChanges×FLAIR
  cross(Original only): 其余 8 个 ROI×序列组合
- 输出长表 radiomics_features_long.csv; 断点续跑按 (patient) 粒度
"""
import os, sys, csv, traceback, logging, gc
import numpy as np
import SimpleITK as sitk
import nibabel as nib

BASE = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts'
INP = r'D:\Raidionics_data\input\2'
REG = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\input_reg\2'   # D盘已满, 配准缓存改落 C 盘
MAIN = r'D:\Raidionics_data\patients\2'
NEC = r'D:\Raidionics_data\patients_nec2\2'
OUT_CSV = os.path.join(BASE, 'radiomics_features_long.csv')
FAIL_CSV = os.path.join(BASE, 'feature_failures.csv')
COHORT = os.environ.get('EF_COHORT', os.path.join(BASE, 'cohort_all_103.csv'))

PRIMARY = [('TumorCore', 'T1CE'), ('Necrosis', 'T1CE'), ('Edema', 'FLAIR'), ('FLAIRChanges', 'FLAIR')]
CROSS = [('TumorCore', 'FLAIR'), ('TumorCore', 'T1w'), ('Necrosis', 'FLAIR'), ('Necrosis', 'T1w'),
         ('Edema', 'T1CE'), ('Edema', 'T1w'), ('FLAIRChanges', 'T1CE'), ('FLAIRChanges', 'T1w')]
SEQ_FILE = {'T1CE': 't1-gd', 'T1w': 't1w', 'FLAIR': 't2-flair'}


def find(d, kw, exclude=None):
    if not d or not os.path.isdir(d):
        return None
    kw = kw.lower()
    for f in sorted(os.listdir(d)):
        fl = f.lower()
        if kw in fl and (exclude is None or exclude.lower() not in fl):
            return os.path.join(d, f)
    return None


def geom(img):
    return (img.GetSize(), img.GetSpacing(), img.GetOrigin(), img.GetDirection())


def geom_same(a, b, tol=1e-4):
    sa, sb = geom(a), geom(b)
    if sa[0] != sb[0]:
        return False
    return all(np.allclose(np.array(x), np.array(y), rtol=0, atol=tol) for x, y in zip(sa[1:], sb[1:]))


def align_geom(mask_img, img):
    """几何仅差机器精度时, 把掩码几何强制对齐到影像几何。
    ITK 的 LabelStatisticsImageFilter 用 1e-6 容差判"同一物理空间",
    nibabel/sitk 浮点往返误差(1e-5级)会让其报错, CopyInformation 可精确消除。"""
    try:
        if mask_img.GetSize() == img.GetSize():
            mask_img.CopyInformation(img)
            return True
    except Exception:
        pass
    return False


def make_mask_image(arr_nib, ref_img):
    """nibabel [x,y,z] uint8 -> sitk image, 几何取自 ref"""
    a = np.asanyarray(arr_nib.dataobj)
    a = np.transpose(a, (2, 1, 0)).astype(np.uint8)
    im = sitk.GetImageFromArray(a)
    im.CopyInformation(ref_img)
    return im


def load_mask(tree, patient, kw, ref_img):
    f = find(os.path.join(tree, patient, 'T0'), kw)
    if f is None:
        return None, None
    return make_mask_image(nib.load(f), ref_img), f


def load_disp(wrp_fn):
    """读取 ANTs warp 场为位移变换。
    统一用 float64 位移场, 与既有 78 例结果保持同一实现(方法学一致)。"""
    img = sitk.ReadImage(wrp_fn)
    disp = sitk.DisplacementFieldTransform(sitk.Cast(img, sitk.sitkVectorFloat64))
    del img
    return disp


def get_sequence(patient, seq, ref_img):
    """返回 (sitk image, 来源标记 raw/registered)"""
    raw_fn = find(os.path.join(INP, patient, 'T0'), SEQ_FILE[seq])
    if raw_fn is None:
        return None, 'missing'
    raw = sitk.ReadImage(raw_fn)
    if seq == 'T1CE':
        return raw, 'raw'
    if geom_same(raw, ref_img):
        return raw, 'raw'
    # 重采样(带缓存)
    cache = os.path.join(REG, patient, SEQ_FILE[seq] + '_reg.nii.gz')
    if os.path.exists(cache):
        return sitk.ReadImage(cache), 'registered'
    trdir = os.path.join(MAIN, patient, 'Transforms')
    tr = None
    for d in sorted(os.listdir(trdir)):
        dl = d.lower()
        if (SEQ_FILE[seq] in dl) and ('t1-gd' in dl) and ('to-' in dl) and ('mni' not in dl):
            tr = os.path.join(trdir, d)
            break
    aff_fn = find(tr, 'genericaffine.mat', exclude='inverse')
    wrp_fn = find(tr, 'warp.nii.gz', exclude='inverse')
    if aff_fn is None or wrp_fn is None:
        return None, 'no_transform'
    aff = sitk.ReadTransform(aff_fn)
    comp = sitk.CompositeTransform([aff, load_disp(wrp_fn)])
    res = sitk.Resample(sitk.Cast(raw, sitk.sitkFloat32), ref_img, comp,
                        sitk.sitkLinear, 0.0, sitk.sitkFloat32)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    sitk.WriteImage(res, cache, useCompression=True)
    return res, 'registered'


def build_extractor(with_wavelet, with_shape):
    from radiomics import featureextractor
    ex = featureextractor.RadiomicsFeatureExtractor(binWidth=25)
    ex.settings['correctMask'] = True      # 兜底: 残余几何差异自动重采样而非丢样本
    ex.settings['geometryTolerance'] = 1e-3
    ex.disableAllImageTypes()
    ex.enableImageTypeByName('Original')
    if with_wavelet:
        ex.enableImageTypeByName('Wavelet')
    ex.disableAllFeatures()
    for cls in (['shape'] if with_shape else []) + ['firstorder', 'glcm', 'glrlm', 'glszm', 'gldm', 'ngtdm']:
        ex.enableFeatureClassByName(cls)
    return ex


def nib_to_sitk(path):
    """nibabel 读取 NIfTI -> 正确几何的 sitk image (RAS->LPS 转换)"""
    nii = nib.load(path)
    a = np.asanyarray(nii.dataobj)
    a = np.transpose(a, (2, 1, 0))
    im = sitk.GetImageFromArray(a.astype(np.uint8) if a.dtype == np.uint8 else a)
    zooms = np.asarray(nii.header.get_zooms()[:3], dtype=float)
    F = np.diag([-1.0, -1.0, 1.0, 1.0])
    alps = F @ nii.affine
    D = alps[:3, :3] @ np.diag(1.0 / zooms)
    im.SetSpacing(tuple(float(z) for z in zooms))
    im.SetOrigin(tuple(float(x) for x in alps[:3, 3]))
    im.SetDirection(tuple(float(x) for x in D.flatten()))
    return im


def get_mask_space_and_image(patient, roi, seq, ref_img, raws, seqs_t1gd):
    """按掩码几何自动匹配影像:
    - 掩码与 t1gd 网格一致 -> 用该序列的 t1gd 网格版本(原始或配准)
    - 掩码与某序列原生网格一致且 seq 相同 -> 直接用原生影像
    - 否则: 用 ANTs 变换把掩码最近邻重采样到 t1gd 网格
    返回 (mask_img, img, src_note) 或 (None, None, reason)
    """
    tree, kw = (NEC, 'annotation-Necrosis') if roi == 'Necrosis' else \
               (MAIN, {'TumorCore': 'annotation-Tumor_MRI_TumorCore',
                       'Edema': 'annotation-Edema',
                       'FLAIRChanges': 'annotation-FLAIRChanges'}[roi])
    mf = find(os.path.join(tree, patient, 'T0'), kw)
    if mf is None:
        return None, None, 'mask missing'
    mimg_native = nib_to_sitk(mf)
    nz = int((sitk.GetArrayFromImage(mimg_native) > 0).sum())
    if nz < 30:
        return None, None, 'mask too small (%d voxels)' % nz

    # 候选: t1gd 网格版本
    img_t1gd_grid, src = seqs_t1gd.get(seq, (None, ''))
    candidates = []
    if img_t1gd_grid is not None:
        candidates.append((img_t1gd_grid, '%s(%s)' % (seq, src)))
    # 掩码与自身原生序列一致(FLAIR/T1w 原生网格)
    if seq in ('FLAIR', 'T1w') and raws.get(seq) is not None:
        candidates.append((raws[seq], '%s(native)' % seq))

    for img, note in candidates:
        if geom_same(mimg_native, img):
            align_geom(mimg_native, img)
            return mimg_native, img, note
    # 掩码重采样到 t1gd 网格(最近邻)
    if img_t1gd_grid is not None:
        trdir = os.path.join(MAIN, patient, 'Transforms')
        tr = None
        # 掩码若在 flair 原生空间, 用 flair->t1gd 变换; 在 t1w 原生空间同理
        mask_space = None
        for s in ('FLAIR', 'T1w'):
            if raws.get(s) is not None and geom_same(mimg_native, raws[s]):
                mask_space = s
                break
        if mask_space is not None:
            kw2 = SEQ_FILE[mask_space]
            for d in sorted(os.listdir(trdir)):
                dl = d.lower()
                if (kw2 in dl) and ('t1-gd' in dl) and ('to-' in dl) and ('mni' not in dl):
                    tr = os.path.join(trdir, d)
                    break
            aff_fn = find(tr, 'genericaffine.mat', exclude='inverse')
            wrp_fn = find(tr, 'warp.nii.gz', exclude='inverse')
            if aff_fn and wrp_fn:
                aff = sitk.ReadTransform(aff_fn)
                comp = sitk.CompositeTransform([aff, load_disp(wrp_fn)])
                mres = sitk.Resample(mimg_native, ref_img, comp, sitk.sitkNearestNeighbor, 0, sitk.sitkUInt8)
                align_geom(mres, img_t1gd_grid)
                return mres, img_t1gd_grid, '%s(mask resampled from %s)' % (seq, mask_space)
        return None, None, 'mask space unresolvable'
    return None, None, 'seq unavailable'


def process_patient(patient):
    rows = []   # (patient, roi, seq, feature, value)
    fails = []  # (patient, roi, seq, error)
    log = []
    try:
        t1gd_fn = find(os.path.join(INP, patient, 'T0'), 't1-gd')
        if t1gd_fn is None:
            return rows, ['%s: no t1-gd' % patient], []
        ref_img = sitk.ReadImage(t1gd_fn)

        raws = {}
        for seq in ('T1w', 'FLAIR'):
            fn = find(os.path.join(INP, patient, 'T0'), SEQ_FILE[seq])
            raws[seq] = sitk.ReadImage(fn) if fn else None

        seqs_t1gd = {}
        for seq in ('T1CE', 'T1w', 'FLAIR'):
            im, src = get_sequence(patient, seq, ref_img)
            seqs_t1gd[seq] = (im, src)
            if im is None:
                log.append('%s: seq %s unavailable (%s)' % (patient, seq, src))

        for roi, seq in PRIMARY + CROSS:
            try:
                mimg, img, note = get_mask_space_and_image(patient, roi, seq, ref_img, raws, seqs_t1gd)
            except Exception as e:
                # 单对失败不应拖垮整例(此前会返回不完整数据), 记为该对失败后继续
                fails.append((patient, roi, seq, 'pair error: ' + str(e)[:130]))
                continue
            if mimg is None:
                if 'too small' not in note:
                    log.append('%s: %s x %s skipped (%s)' % (patient, roi, seq, note))
                else:
                    fails.append((patient, roi, seq, note))
                continue
            is_primary = (roi, seq) in PRIMARY
            try:
                with_shape = is_primary and seq in ('T1CE', 'FLAIR')
                ex = build_extractor(is_primary, with_shape)
                result = ex.execute(img, mimg, label=1)
                for k, v in result.items():
                    if k.startswith('diagnostics'):
                        continue
                    try:
                        val = float(v)
                    except Exception:
                        continue
                    if not np.isfinite(val):
                        continue
                    rows.append((patient, roi, seq, k, val))
            except Exception as e:
                fails.append((patient, roi, seq, str(e)[:150]))
        # 释放大体积影像, 降低内存碎片(大网格病例峰值可达数 GB)
        for k in list(seqs_t1gd):
            seqs_t1gd[k] = (None, '')
        for k in list(raws):
            raws[k] = None
        ref_img = None
        gc.collect()
        return rows, log, fails
    except Exception:
        log.append('%s: FATAL %s' % (patient, traceback.format_exc()[-300:]))
        return rows, log, fails


def main():
    logging.getLogger('radiomics').setLevel(logging.ERROR)
    cohort = [r['patient'] for r in csv.DictReader(open(COHORT, encoding='utf-8-sig'))]
    done = set()
    if os.path.exists(OUT_CSV):
        with open(OUT_CSV, encoding='utf-8-sig') as f:
            for r in csv.reader(f):
                if r and r[0] != 'patient':
                    done.add(r[0])
    todo = [p for p in cohort if p not in done]
    print('cohort=%d done=%d todo=%d' % (len(cohort), len(done), len(todo)), flush=True)
    if not todo:
        print('ALL_DONE', flush=True)
        return

    new_file = (not os.path.exists(OUT_CSV)) or os.path.getsize(OUT_CSV) == 0
    fout = open(OUT_CSV, 'a', newline='', encoding='utf-8-sig')
    wout = csv.writer(fout)
    if new_file:
        wout.writerow(['patient', 'roi', 'seq', 'feature', 'value'])
    ffail = open(FAIL_CSV, 'a', newline='', encoding='utf-8-sig')
    wfail = csv.writer(ffail)
    if os.path.getsize(FAIL_CSV) == 0:
        wfail.writerow(['patient', 'roi', 'seq', 'error'])

    from multiprocessing import Pool
    with Pool(4) as pool:
        for i, (rows, log, fails) in enumerate(pool.imap_unordered(process_patient, todo)):
            wout.writerows(rows)
            wfail.writerows(fails)
            fout.flush()
            ffail.flush()
            for m in log:
                print('  ' + m, flush=True)
            print('[%d/%d] %s -> %d features' % (i + 1, len(todo), rows[0][0] if rows else '?', len(rows)), flush=True)
    fout.close()
    ffail.close()
    print('ALL_DONE', flush=True)


if __name__ == '__main__':
    main()
