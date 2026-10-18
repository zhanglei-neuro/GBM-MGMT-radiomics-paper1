# -*- coding: utf-8 -*-
r"""
check_phi.py — 上传前的患者隐私自检（务必在每次 git push 之前运行）

用法:
    python tools/check_phi.py                 # 扫描仓库根目录
    python tools/check_phi.py --names names.txt   # 额外用姓名清单匹配

检测项:
  1. 原始 DICOM 路径片段 (NeuroImagepre / NeuroImagepost / Raidionics_data) —— 提示级
  2. DICOM 目录名中出现中文人名 (ImageFileName_<中文>) —— 严重
  3. 形如 disc_<小写拼音> / NeuroImagepre\<大写拼音> 的病例目录名 —— 严重
  4. 姓名清单逐词匹配 —— 严重（清单默认自动读取 _internal/patient_ids_local.csv
     与 _internal/patient_alias_local.csv，可用 --names 覆盖；找不到时降级为不做匹配）
  5. 疑似患者数据文件 (*.dcm, *.nii*, GBM_clinical_data.xlsx, patient_ids_local.csv 等) —— 严重

注意：`_internal/` 本身存放姓名映射，扫描整包（root='.'）时必然报严重 —— 验收请分别对
`github` 与 `zenodo` 两个上传目录运行，二者均需为「严重 0 项」。

退出码 0 = 未发现严重问题；1 = 发现严重问题（不要 push）。
"""
import os, re, sys, argparse

SEVERE, WARN = [], []

CJK_NAME = re.compile(r'ImageFileName[_\s]*[\u4e00-\u9fff][\u4e00-\u9fff\s]*')
DICOM_DIR = re.compile(r'(?:NeuroImagepre|NeuroImagepost|Raidionics_data)[\\/]+'
                       r'([A-Za-z][A-Za-z]{3,})')
DISC_DIR = re.compile(r'\bdisc[_-]([a-z][a-z]{3,})\b')
DICOM_ANY = re.compile(r'(NeuroImagepre|NeuroImagepost|Raidionics_data)', re.I)
BAD_EXT = ('.dcm', '.nii', '.nii.gz', '.ima', '.dicom')
BAD_NAME = ('GBM_clinical_data.xlsx', 'cohort_MGMT.csv', 'id_mapping',
            'series_selection.csv',
            # 患者编号 ↔ 姓名对照表（含 PHI）。正常情况下放在 _internal/，
            # 但复制到仓库内即可导致泄露 —— 必须在此拦截。
            'patient_ids_local.csv', 'patient_alias_local.csv')
GENERIC = {'data', 'input', 'output', 'patients', 'patients_nec2', 'disc',
           'study', 'series', 'image', 'patient', 'pre', 'post', 't0', 't1',
           # 变量名/字段名，非病例目录名（'disc_desc' 会命中 DISC_DIR）
           'desc', 'desc_map', 'seriesdescription'}

TEXT_EXT = ('.py', '.md', '.txt', '.json', '.csv', '.yml', '.yaml', '.cfg',
            '.ini', '.rst')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root', nargs='?', default='.')
    ap.add_argument('--names', help='姓名清单文件（每行一个）')
    a = ap.parse_args()

    names = []
    if a.names and os.path.isfile(a.names):
        names = [l.strip() for l in open(a.names, encoding='utf-8') if len(l.strip()) >= 4]
    else:
        # 未显式指定清单时，自动读取本地映射表；否则本工具无法做姓名匹配（历史盲区）。
        _here = os.path.dirname(os.path.abspath(__file__))
        for cand, col in ((os.path.join(_here, '..', '..', '_internal', 'patient_ids_local.csv'), 1),
                          (os.path.join(_here, 'patient_ids_local.csv'), 1),
                          (os.path.join(_here, '..', '..', '_internal', 'patient_alias_local.csv'), 0),
                          (os.path.join(_here, 'patient_alias_local.csv'), 0)):
            if os.path.isfile(cand):
                rows = [l.strip().split(',') for l in
                        open(cand, encoding='utf-8-sig') if l.strip()]
                names += [r[col].strip() for r in rows[1:]
                          if len(r) > col and len(r[col].strip()) >= 4]
        if names:
            print('[info] 姓名清单自动加载 %d 个（含拼写变体）' % len(set(names)))
        else:
            print('[warn] 未找到姓名清单（--names / _internal/patient_ids_local.csv），'
                  '本次不做姓名匹配')
    name_pats = [(n, re.compile(re.escape(n), re.I)) for n in
                 sorted(names, key=len, reverse=True)]

    for dp, dirs, fs in os.walk(a.root):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', '.venv')]
        for fn in fs:
            p = os.path.join(dp, fn)
            rel = os.path.relpath(p, a.root).replace('\\', '/')

            if fn.lower().endswith(BAD_EXT) or (
                    any(b in fn for b in BAD_NAME)
                    and not fn.lower().endswith('.example')):
                SEVERE.append((rel, '疑似患者数据文件'))
                continue

            if not fn.lower().endswith(TEXT_EXT):
                continue
            try:
                s = open(p, encoding='utf-8', errors='ignore').read()
            except Exception:
                continue

            if CJK_NAME.search(s):
                SEVERE.append((rel, '中文人名出现在 DICOM 目录名'))
            for m in set(DICOM_DIR.findall(s)) | set(DISC_DIR.findall(s)):
                if m.lower() not in GENERIC and not m.startswith('<'):
                    SEVERE.append((rel, '疑似病例目录名：%s' % m))
            for n, pt in name_pats:
                if pt.search(s):
                    SEVERE.append((rel, '姓名清单命中：%s' % n))
            if DICOM_ANY.search(s):
                WARN.append((rel, '含原始 DICOM 路径（无姓名则为提示级）'))

    print('=' * 66)
    print('患者隐私自检报告')
    print('=' * 66)
    if not SEVERE and not WARN:
        print('未发现问题，可以 push。')
    for rel, why in sorted(set(SEVERE)):
        print('  [严重] %-58s %s' % (rel, why))
    for rel, why in sorted(set(WARN)):
        print('  [提示] %-58s %s' % (rel, why))
    print('-' * 66)
    print('严重 %d 项 / 提示 %d 项' % (len(set(SEVERE)), len(set(WARN))))
    if SEVERE:
        print('\n发现严重问题：请先运行 tools/sanitize_phi.py 或手工处理后再 push。')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
