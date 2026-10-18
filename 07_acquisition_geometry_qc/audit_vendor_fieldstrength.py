# -*- coding: utf-8 -*-
r"""全队列场强/厂商核验 —— 逐例读取 E:\NeuroImagepre\<患者>\ 下的首个可读 DICOM 头，
输出 Manufacturer / ManufacturerModelName / MagneticFieldStrength。

背景：临床表 GBM_clinical_data.xlsx 的 MRI_vendor 列存在陈旧/错误标注
（如 P055 被标为 UIH，实为 Siemens），不能作为平台归属依据；
Paper 1 已改用 DICOM 实测厂商。本脚本将该做法推广到全部患者，
并检出 Paper 2 随访队列中混入的第 2 例 United Imaging uMR560 1.5 T（P055）。

用法: python audit_vendor_fieldstrength.py <患者名单csv> <输出csv>
      患者名单 csv 需含列 'patient'
"""
import os
import sys
import pandas as pd
import pydicom

DICOM_ROOT = r'E:\NeuroImagepre'


def scan_first_hdr(base):
    """深度遍历目录，返回首个可读 DICOM 头的 (Manufacturer, Model, FieldStrength)。"""
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
        print('\n非 3.0 T 病例（应排除）:')
        print(bad.to_string())


if __name__ == '__main__':
    main()
