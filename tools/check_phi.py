# -*- coding: utf-8 -*-
r"""
check_phi.py -- patient privacy self-check before upload (run this before every git push)

Usage:
    python tools/check_phi.py                     # scan the repository root
    python tools/check_phi.py --names names.txt   # additionally match against a name list

What it detects:
  1. Fragments of original DICOM paths (NeuroImagepre / NeuroImagepost / Raidionics_data) -- notice level
  2. Chinese personal names inside DICOM directory names (ImageFileName_<Chinese>) -- severe
  3. Case directory names of the form disc_<lower-case pinyin> / NeuroImagepre\<upper-case pinyin> -- severe
  4. Word-by-word matching against the name list -- severe (the list is picked up automatically from
     _internal/patient_ids_local.csv and _internal/patient_alias_local.csv; --names overrides it, and
     matching degrades to a no-op if none is found)
  5. Files that look like patient data (*.dcm, *.nii*, GBM_clinical_data.xlsx, patient_ids_local.csv, ...) -- severe

Note: `_internal/` itself holds the name mapping, so scanning the whole package (root='.') will always report severe
findings. For sign-off, run it against the two upload directories, `github` and `zenodo`, separately; both must
report "0 severe".

Exit code 0 = no severe findings; 1 = severe findings (do not push).
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
            # patient code <-> name mapping table (contains PHI). It normally lives in _internal/,
            # but copying it into the repository is enough to leak it, so it must be blocked here.
            'patient_ids_local.csv', 'patient_alias_local.csv')
GENERIC = {'data', 'input', 'output', 'patients', 'patients_nec2', 'disc',
           'study', 'series', 'image', 'patient', 'pre', 'post', 't0', 't1',
           # variable and field names, not case directory names ('disc_desc' would otherwise hit DISC_DIR)
           'desc', 'desc_map', 'seriesdescription'}

TEXT_EXT = ('.py', '.md', '.txt', '.json', '.csv', '.yml', '.yaml', '.cfg',
            '.ini', '.rst')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('root', nargs='?', default='.')
    ap.add_argument('--names', help='name list file (one name per line)')
    a = ap.parse_args()

    names = []
    if a.names and os.path.isfile(a.names):
        names = [l.strip() for l in open(a.names, encoding='utf-8') if len(l.strip()) >= 4]
    else:
        # without an explicit list, read the local mapping table; otherwise this tool cannot match names at all (a historical blind spot)
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
            print('[info] name list auto-loaded, %d entries (including spelling variants)' % len(set(names)))
        else:
            print('[warn] no name list found (--names / _internal/patient_ids_local.csv), '
                  'so no name matching will be done this run')
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
                SEVERE.append((rel, 'file that may contain patient data'))
                continue

            if not fn.lower().endswith(TEXT_EXT):
                continue
            try:
                s = open(p, encoding='utf-8', errors='ignore').read()
            except Exception:
                continue

            if CJK_NAME.search(s):
                SEVERE.append((rel, 'Chinese personal name in a DICOM directory name'))
            for m in set(DICOM_DIR.findall(s)) | set(DISC_DIR.findall(s)):
                if m.lower() not in GENERIC and not m.startswith('<'):
                    SEVERE.append((rel, 'possible case directory name: %s' % m))
            for n, pt in name_pats:
                if pt.search(s):
                    SEVERE.append((rel, 'matched the name list: %s' % n))
            if DICOM_ANY.search(s):
                WARN.append((rel, 'contains an original DICOM path (notice level unless a name is present)'))

    print('=' * 66)
    print('Patient privacy self-check report')
    print('=' * 66)
    if not SEVERE and not WARN:
        print('No problems found; safe to push.')
    for rel, why in sorted(set(SEVERE)):
        print('  [severe] %-58s %s' % (rel, why))
    for rel, why in sorted(set(WARN)):
        print('  [notice] %-58s %s' % (rel, why))
    print('-' * 66)
    print('%d severe / %d notice' % (len(set(SEVERE)), len(set(WARN))))
    if SEVERE:
        print('\nSevere problems found: run tools/sanitize_phi.py or fix them by hand before pushing.')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
