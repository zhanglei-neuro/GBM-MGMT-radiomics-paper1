# -*- coding: utf-8 -*-
r"""步骤1b: 从 dicom_inventory.csv 重新做角色判定与序列选择 (修正版)
修正点:
  - FLAIR 识别: 'flair' / 'dark-fluid' / 'dark_fluid' / ('t2'+'tirm')  (西门子协议 dark-fluid 即 FLAIR)
  - 排除干扰序列: localizer/scout/灌注(perf)/弥散(adc,diff,trace)/MPR重建/Range
  - T1-CE 排序: 无缺层 > 层数多(3D优于2D) > 轴位(tra)
产出: series_selection.csv (每个患者选定的三条序列目录)
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'tools'))
from patient_ids import to_name   # noqa: E402  患者编号 <-> 本地姓名
import csv, glob, os, re

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dicom_inventory.csv')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'series_selection.csv')

EXCLUDE_TOKENS = ['localizer', 'scout', 'perf', 'adc', 'diff', 'trace', 'range', 'mpr 1mm']

def dl(s):
    return s.lower()

def excluded(desc):
    d = dl(desc)
    # GE 序列名 t1_mpr_sag_fs_iso+C 是 3D T1-CE 本体, 不算 MPR 重建
    is_ge_seq = d.startswith('t1_mpr_')
    if not is_ge_seq and (re.search(r'mpr[_ ](cor|tra|sag)', d) or '<mpr' in d):
        return True
    for t in EXCLUDE_TOKENS:
        if t in d:
            return True
    return False

def is_flair_desc(desc):
    d = dl(desc)
    if 'flair' in d:
        return True
    if ('dark-fluid' in d) or ('dark_fluid' in d) or ('darkfluid' in d):
        return True
    if ('t2' in d) and ('tirm' in d):
        return True
    return False

def classify(r):
    d = dl(r['desc'])
    agent = bool(r['contrast_agent'].strip())
    if excluded(r['desc']):
        return 'other'
    has_t1 = ('t1' in d) or ('mprage' in d)
    if has_t1 and agent:
        # 铁律: ContrastBolusAgent 非空优先于名字 (GE 的 t1_fse_flair_+c 是增强 T1)
        return 't1ce'
    if is_flair_desc(r['desc']) and not agent and not ('t1' in d):
        return 'flair'
    if has_t1 and not agent:
        return 't1'
    return 'other'

def rank(r):
    tra = 0 if 'tra' in dl(r['desc']) else 1
    return (r['gap_bad'] == 'True', -int(r['n_files']), tra)

def main():
    rows = list(csv.DictReader(open(SRC, encoding='utf-8-sig')))
    patients = sorted(set(r['patient'] for r in rows))
    # 定向覆盖: t1_mpr_sag_fs_iso+C 混入两个错位体积副本的患者, 改用干净 _ND 变体
    # 键为去标识编号; 患者目录名以姓名命名, 故路径经 to_name() 还原为本地目录名
    PRE_ROOT = os.environ.get('PATIENT_ROOT', r'D:\NeuroImagepre')

    def _ovr_dir(pid, date, ser):
        base = os.path.join(PRE_ROOT, to_name(pid), date)
        sub = sorted(glob.glob(os.path.join(base, 'ImageFileName_*')))
        return os.path.join(sub[0] if sub else base, ser)

    OVERRIDES = {pid: ('t1ce', _ovr_dir(pid, date, ser)) for pid, date, ser in
                 [('P032', '20231228', '8'),
                  ('P051', '20220314', '8'),
                  ('P108', '20230506', '12')]}
    sel_rows = []
    for p in patients:
        prs = [r for r in rows if r['patient'] == p]
        for r in prs:
            r['role'] = classify(r)
        sel = {'patient': p}
        for role in ('t1ce', 'flair', 't1'):
            cand = [r for r in prs if r['role'] == role]
            cand.sort(key=rank)
            c = cand[0] if cand else None
            if p in OVERRIDES and OVERRIDES[p][0] == role:
                ovr_dir = OVERRIDES[p][1]
                ovr = [r for r in prs if r['series_dir'] == ovr_dir]
                if ovr:
                    c = ovr[0]
                    sel[role + '_desc'] = f"OVERRIDE {c['desc']} ({c['n_files']}slices sp={c['slice_sp']})"
            else:
                sel[role + '_desc'] = (f"{c['desc']} ({c['n_files']}slices sp={c['slice_sp']})") if c else ''
            sel[role] = c['series_dir'] if c else ''
            sel[role + '_gap'] = ('BAD:' + c['gap_report']) if (c and c['gap_bad'] == 'True') else ('ok' if c else '-')
        sel_rows.append(sel)
    with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['patient', 't1ce', 't1ce_desc', 't1ce_gap',
                                          'flair', 'flair_desc', 'flair_gap',
                                          't1', 't1_desc', 't1_gap'])
        w.writeheader()
        w.writerows(sel_rows)
    n_ok = sum(1 for s in sel_rows if s['t1ce'] and s['flair'] and s['t1'])
    print(f"total={len(sel_rows)} complete={n_ok}")
    for s in sel_rows:
        flags = ''.join(['1' if s['t1ce'] else '0', '1' if s['flair'] else '0', '1' if s['t1'] else '0'])
        mark = 'BADGAP' if (s['t1ce_gap'].startswith('BAD') or s['flair_gap'].startswith('BAD') or s['t1_gap'].startswith('BAD')) else ''
        print(f"{s['patient']}: {flags} {mark} | CE={s['t1ce_desc']} | FL={s['flair_desc']} | T1={s['t1_desc']}")

if __name__ == '__main__':
    main()
