# -*- coding: utf-8 -*-
r"""configure_paths.py -- rewrite the original absolute paths in the scripts to point at your own directories

Usage (run --dry first to preview):
    python tools/configure_paths.py --dry
    python tools/configure_paths.py --set "D:\\my\\clinical.xlsx=D:\\data\\clinical.xlsx"

Without --set it goes into interactive mode and asks about each path in turn.
Replacement only touches .py files and matches both the forward-slash and backslash spellings.
"""
import os, re, sys, argparse, collections

DEFAULTS = [
    r'E:\NeuroImagepre',
    r'D:\NeuroImagepre',
    r'D:\NeuroImagepost',
    r'D:\Raidionics_data',
    r'D:\1桌面文件\AI\GBM_clinical_data.xlsx',
    r'D:\1桌面文件\AI\Glioma_Patients_Classified_merged.xlsx',
    r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57',
    r'D:\WorkBuddy\2026-10-06-10-35-00\outputs\paper1',
    r'D:\Anaconda\envs\nnunet',
]


def apply_rules(root, rules, dry=True):
    n_file = n_hit = 0
    for dp, dirs, fs in os.walk(root):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', '.venv')]
        for fn in fs:
            if not fn.endswith('.py'):
                continue
            p = os.path.join(dp, fn)
            s0 = open(p, encoding='utf-8', errors='ignore').read()
            s, hits = s0, 0
            for old, new in rules:
                for variant in {old, old.replace('\\', '/'), old.replace('/', '\\')}:
                    if variant in s:
                        s = s.replace(variant, new)
                        hits += s0.count(variant)
            if hits and s != s0:
                n_file += 1
                n_hit += 1
                print('  %-58s rewritten' % os.path.relpath(p, root))
                if not dry:
                    open(p, 'w', encoding='utf-8').write(s)
    print('\n%s: %d files' % ('dry run' if dry else 'done', n_file))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--set', action='append', default=[],
                    metavar='OLD=NEW')
    ap.add_argument('--dry', action='store_true')
    a = ap.parse_args()

    if a.set:
        rules = []
        for kv in a.set:
            if '=' not in kv:
                print('ignored (no = sign):', kv)
                continue
            old, new = kv.split('=', 1)
            rules.append((old.strip(), new.strip()))
    else:
        print('no --set given, entering interactive mode. Press Enter to skip an entry.\n')
        rules = []
        for old in DEFAULTS:
            new = input('  %s\n    → ' % old).strip()
            if new:
                rules.append((old, new))
        if not rules:
            print('no paths configured, exiting.')
            return

    print('\nthe following replacements will be made:')
    for old, new in rules:
        print('  %s  →  %s' % (old, new))
    print()
    apply_rules(a.root, rules, dry=a.dry)
    if a.dry:
        print('hint: drop --dry to actually write.')


if __name__ == '__main__':
    main()
