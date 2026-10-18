# -*- coding: utf-8 -*-
r"""
sanitize_phi.py — 把脚本中的患者姓名替换为 <PATIENT> 占位符

用法:
    python tools/sanitize_phi.py --names names.txt            # 用清单匹配（推荐）
    python tools/sanitize_phi.py --names names.txt --root .
    python tools/sanitize_phi.py --names names.txt --dry       # 只报告不修改

清单文件格式：每行一个姓名（拼音或中文），例如
    <Surname><Given>
    <Surname><Given>
会自动忽略长度 < 4 的条目。替换时按长度降序处理，避免短名误伤长名。
中文姓名在 `ImageFileName_<中文>` 上下文中会被整体替换。
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
    print('姓名清单 %d 条' % len(pats))

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
                print('  %-60s %3d 处' % (os.path.relpath(p, a.root), hits))
                if not a.dry:
                    open(p, 'w', encoding='utf-8').write(s)
    print('\n%s：%d 个文件 / %d 处替换'
          % ('预演' if a.dry else '完成', total_f, total_h))


if __name__ == '__main__':
    sys.exit(main())
