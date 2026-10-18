# -*- coding: utf-8 -*-
r"""
configure_paths.py — 把脚本中的原始绝对路径改写为你自己的目录

用法（推荐先 --dry 预演）:
    python tools/configure_paths.py --dry
    python tools/configure_paths.py \
        --set "E:\NeuroImagepre=D:\data\NeuroImagepre" \
        --set "D:\Raidionics_data=E:\Raidionics\out" \
        --set "D:\1桌面文件\AI\GBM_clinical_data.xlsx=D:\data\clinical.xlsx"

不带 --set 时进入交互模式，逐条询问是否替换。
替换仅作用于 .py 文件，且同时匹配正斜杠 / 反斜杠两种写法。
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
                print('  %-58s 已改写' % os.path.relpath(p, root))
                if not dry:
                    open(p, 'w', encoding='utf-8').write(s)
    print('\n%s：%d 个文件' % ('预演' if dry else '完成', n_file))


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
                print('忽略（缺少 =）：', kv)
                continue
            old, new = kv.split('=', 1)
            rules.append((old.strip(), new.strip()))
    else:
        print('未提供 --set，进入交互模式。直接回车 = 跳过该条。\n')
        rules = []
        for old in DEFAULTS:
            new = input('  %s\n    → ' % old).strip()
            if new:
                rules.append((old, new))
        if not rules:
            print('未配置任何路径，退出。')
            return

    print('\n将执行以下替换：')
    for old, new in rules:
        print('  %s  →  %s' % (old, new))
    print()
    apply_rules(a.root, rules, dry=a.dry)
    if a.dry:
        print('提示：去掉 --dry 即真正写入。')


if __name__ == '__main__':
    main()
