# -*- coding: utf-8 -*-
"""阶段0收尾: 锁定 MGMT 建模队列
- 以 segmentation_timing.csv status=OK 的 103 例为基准
- 关联合并表(标签/厂商/临床), 校验掩码齐全度(TumorCore/Edema/FLAIRChanges + Necrosis新掩码)
- 输出 cohort_all_103.csv 与 cohort_MGMT.csv
"""
import os, csv
import pandas as pd
import SimpleITK as sitk

BASE = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\scripts'
MAIN = r'D:\Raidionics_data\patients\2'
NEC = r'D:\Raidionics_data\patients_nec2\2'
MERGED = r'D:\1桌面文件\AI\Glioma_Patients_Classified_merged.xlsx'
OUT_ALL = os.path.join(BASE, 'cohort_all_103.csv')
OUT_MGMT = os.path.join(BASE, 'cohort_MGMT.csv')

timing = list(csv.DictReader(open(os.path.join(BASE, 'segmentation_timing.csv'), encoding='utf-8-sig')))
ok_patients = [r['patient'] for r in timing if r['status'] == 'OK']

df = pd.read_excel(MERGED, sheet_name='GBM_merged')
df['Name'] = df['Name'].astype(str).str.strip()
mdf = df.set_index('Name')

def find_mask(p, tree, patterns):
    t0 = os.path.join(tree, p, 'T0')
    if not os.path.isdir(t0):
        return None
    files = os.listdir(t0)
    for pat in patterns:
        for f in files:
            if pat in f and f.endswith('.nii.gz'):
                return os.path.join(t0, f)
    return None

TC_PAT = ['annotation-Tumor_MRI_TumorCore']
ED_PAT = ['annotation-Edema']
FL_PAT = ['annotation-FLAIRChanges']
NEC_PAT = ['annotation-Necrosis']

rows = []
for p in ok_patients:
    r = {'patient': p}
    if p in mdf.index:
        m = mdf.loc[p]
        if hasattr(m, 'iloc') and getattr(m, 'ndim', 1) > 1:
            m = m.iloc[0]
        r['MGMT'] = str(m.get('MGMT methylation', '')).strip()
        r['IDH'] = str(m.get('IDH1/2', '')).strip()
        r['TERT'] = str(m.get('TERT', '')).strip()
        r['vendor'] = str(m.get('MRI_vendor', '')).strip()
        r['age'] = m.get('Age', '')
        r['gender'] = m.get('Gender', '')
        r['site'] = m.get('Site', '')
        r['surgery_date'] = str(m.get('Surgery date', '')).strip()
    else:
        r.update({'MGMT': '', 'IDH': '', 'TERT': '', 'vendor': '', 'age': '', 'gender': '',
                  'site': '', 'surgery_date': ''})
    masks = {
        'TumorCore': find_mask(p, MAIN, TC_PAT),
        'Edema': find_mask(p, MAIN, ED_PAT),
        'FLAIRChanges': find_mask(p, MAIN, FL_PAT),
        'Necrosis': find_mask(p, NEC, NEC_PAT),
    }
    r['mask_ok'] = all(v for v in masks.values())
    r['missing_mask'] = ','.join(k for k, v in masks.items() if not v)
    # 空坏死掩码记录(体积=0 属合法结局, 但单独标记)
    if masks['Necrosis']:
        a = sitk.GetArrayFromImage(sitk.ReadImage(masks['Necrosis']))
        r['necrosis_voxel'] = int((a > 0).sum())
    else:
        r['necrosis_voxel'] = ''
    rows.append(r)

fields = ['patient', 'MGMT', 'IDH', 'TERT', 'vendor', 'age', 'gender', 'site',
          'surgery_date', 'mask_ok', 'missing_mask', 'necrosis_voxel']
with open(OUT_ALL, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

mgmt = [r for r in rows if r['MGMT'] in ('+', '-', '+1', '-1', '阳性', '阴性')]
with open(OUT_MGMT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(mgmt)

pos = sum(1 for r in mgmt if r['MGMT'].startswith('+'))
neg = sum(1 for r in mgmt if r['MGMT'].startswith('-'))
ven = {}
for r in mgmt:
    ven[r['vendor']] = ven.get(r['vendor'], 0) + 1
mask_bad = [r['patient'] for r in mgmt if not r['mask_ok']]

out = [
    'cohort_all=%d (mask_ok=%d)' % (len(rows), sum(1 for r in rows if r['mask_ok'])),
    'MGMT labeled=%d (+=%d, -=%d)' % (len(mgmt), pos, neg),
    'vendor dist=%s' % ven,
    'mask_missing_in_MGMT_cohort=%s' % mask_bad,
    'files: %s , %s' % (OUT_ALL, OUT_MGMT),
]
open(r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57\.workbuddy\cohort_res.txt', 'w', encoding='utf-8').write('\n'.join(out))
print('done')
