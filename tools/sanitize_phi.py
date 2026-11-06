# -*- coding: utf-8 -*-
r"""sanitize_phi.py -- replace patient names in the scripts with the <PATIENT> placeholder

Usage:
    python tools/sanitize_phi.py --names names.txt            # match against a list (recommended)
    python tools/sanitize_phi.py --names names.txt --dry      # report only, change nothing

List format: one name per line (pinyin or Chinese), for example those exported from the local mapping table.
Entries shorter than 4 characters are ignored. Replacement runs longest-first so that a short name cannot damage a
longer one. A Chinese name inside `ImageFileName_<Chinese>` is replaced as a whole.
"""
import os, re, sys, argparse

CJK_IN_FILE = re.compile(r'(ImageFileName[_\s]*)[\u4e00-\u9fff][\u4e00-\u9fff\s]*')
TOKEN = '<PATIENT>'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--names', required=True)
    ap.add_argument('--root', default='.')
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--ext', default='.py')
    a = ap.parse_args()

    names = [l.strip() for l in open(a.names, encoding='utf-8') if len(l.strip()) >= 4]
    pats = [(n, re.compile(re.escape(n), re.I))
            for n in sorted(set(names), key=len, reverse=True)]
    print('%d names in the list' % len(pats))

    total_f, total_h = 0, 0
    for dp, dirs, fs in os.walk(a.root):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', '.venv')]
        for fn in fs:
            if not fn.endswith(tuple(a.ext.split(','))):
                continue
            p = os.path.join(dp, fn)
            s0 = open(p, encoding='utf-8', errors='ignore').read()
            s, hits = s0, 0
            for n, pt in pats:
                s, k = pt.subn(TOKEN, s)
                hits += k
            s, k = CJK_IN_FILE.subn(r'\1' + TOKEN, s)
            hits += k
            s = re.sub(TOKEN + '+', TOKEN, s)
            if hits:
                total_f += 1
                total_h += hits
                print('  %-60s %3d hits' % (os.path.relpath(p, a.root), hits))
                if not a.dry:
                    open(p, 'w', encoding='utf-8').write(s)
    print('\n%s: %d files / %d replacements'
          % ('dry run' if a.dry else 'done', total_f, total_h))


if __name__ == '__main__':
    sys.exit(main())
