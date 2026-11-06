# -*- coding: utf-8 -*-
r"""Infer the MRI scanner vendor from the series description and cross-validate it against the post-operative ground truth.

Basis for the mapping: the correspondence between the true Manufacturer (0008,0070) in the post-operative
(NeuroImagepost) data and SeriesDescription.
"""
import csv, os, sys, io, json, glob
from collections import Counter, defaultdict
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
BASE = r'C:\Users\zhanglei\WorkBuddy\2026-09-18-22-15-57'

SIEMENS, GE, UIH, OTHER = 'Siemens', 'GE', 'UIH', 'Unknown'

RULES = [
    # ---- Siemens (most specific, so it must be tested first; includes the Siemens perfusion parameter maps ep2d_perf_*_TTP) ----
    ('tirm', SIEMENS), ('qtse', SIEMENS), ('ep2d', SIEMENS), ('tse', SIEMENS),
    ('mprage_sag_p2', SIEMENS), ('aahead', SIEMENS), ('resolve_', SIEMENS),
    ('mpr_sag_fs_iso', SIEMENS), ('dark-fluid', SIEMENS), ('fl2d', SIEMENS),
    ('localizer', SIEMENS), ('mpr range', SIEMENS), ('mpr 1mm range', SIEMENS),
    ('mpr cor 1mm range', SIEMENS), ('mpr thick range', SIEMENS),
    ('abi1_t1iso_mprage', SIEMENS),
    # ---- UIH, United Imaging ----
    ('epi_dwi_tra', UIH), ('gre_scout', UIH), ('fse_flair', UIH), ('t2_fse_tra', UIH),
    ('fse_flair_tra_fs', UIH), ('t1_gre_', UIH),
    # ---- GE ----
    ('fspgr', GE), ('mp-rage', GE), ('propeller', GE), ('muse', GE), ('ssfse', GE),
    ('adc (10^-6', GE), ('5mm ob_', GE), ('probe sv', GE), ('pwi', GE),
    ('rcbv', GE), ('rcbf', GE), ('rmtt', GE), ('mtt', GE), ('ttp', GE),
    ('cbv', GE), ('cbf', GE),
    ('ax t2 fs flair', GE), ('sag t1 flair', GE), ('ax t1 flair', GE),
    ('3-pl loc gre', GE), ('processed images', GE),
]

def vendor_of(desc):
    d = desc.lower()
    for pat, v in RULES:
        if pat in d:
            return v
    return OTHER

def build_preop():
    rows = list(csv.DictReader(open(os.path.join(BASE, 'scripts', 'dicom_inventory.csv'), encoding='utf-8-sig')))
    per = defaultdict(list)
    for r in rows:
        per[r['patient']].append(r['desc'])
    res = {}
    for pat in sorted(per):
        votes = defaultdict(list)
        for d in per[pat]:
            votes[vendor_of(d)].append(d)
        known = {k: sorted(set(v)) for k, v in votes.items() if k != OTHER}
        if known:
            main_v = max(known.items(), key=lambda kv: len(kv[1]))[0]
        else:
            main_v = OTHER
        res[pat] = dict(
            vendor=main_v, n_series=len(per[pat]),
            evidence=known.get(main_v, []),
            minor={k: v for k, v in known.items() if k != main_v},
            unknown=sorted(set(votes.get(OTHER, []))))
    return res

def build_postop_truth():
    import pydicom
    root = r'D:\1桌面文件\AI\NeuroImagepost'
    truth = {}
    for pat in sorted(os.listdir(root)):
        pdir = os.path.join(root, pat)
        if not os.path.isdir(pdir):
            continue
        c = Counter()
        for dirpath, dirnames, filenames in os.walk(pdir):
            dcms = [f for f in filenames if f.lower().endswith('.dcm')]
            if not dcms:
                continue
            try:
                d = pydicom.dcmread(os.path.join(dirpath, dcms[0]), stop_before_pixels=True, force=True)
            except Exception:
                continue
            man = str(getattr(d, 'Manufacturer', '') or '')
            if man:
                c[man] += 1
        if c:
            truth[pat] = c.most_common(1)[0][0]
    return truth

def norm(v):
    s = (v or '').lower()
    if 'siemens' in s:
        return 'Siemens'
    if 'ge ' in s or s.startswith('ge') or 'gemedical' in s.replace(' ', ''):
        return 'GE'
    if 'united' in s or s.strip() == 'uih' or '联影' in s:
        return 'UIH'
    return v or 'Unknown'

def main():
    pre = build_preop()
    post = build_postop_truth()
    out = []
    out.append('== vendor inferred for this pre-operative batch ==  %d cases in total' % len(pre))
    for pat in sorted(pre):
        r = pre[pat]
        flags = []
        if r['minor']:
            flags.append('second choice:' + ','.join('%s(%d)' % (k, len(v)) for k, v in r['minor'].items()))
        if r['unknown']:
            flags.append('unrecognised:' + ','.join(r['unknown']))
        out.append('%-16s %-9s series=%-3d main evidence=%-3d %s' % (
            pat, r['vendor'], r['n_series'], len(r['evidence']), '; '.join(flags)))
    out.append('')
    out.append('pre-operative distribution: %r' % dict(Counter(r['vendor'] for r in pre.values()).most_common()))
    out.append('')
    both = sorted(set(pre) & set(post))
    out.append('== cross-validation: pre-operative inference (naming) vs post-operative truth ==  %d cases in common' % len(both))
    agree = dis = 0
    for pat in both:
        a, b = pre[pat]['vendor'], norm(post[pat])
        ok = (a == b)
        agree += ok
        dis += (not ok)
        out.append('%-16s inferred=%-9s truth=%-9s %s' % (pat, a, b, 'OK' if ok else 'MISMATCH'))
    out.append('')
    out.append('agree %d / disagree %d / agreement %.1f%%' % (agree, dis, 100.0 * agree / max(len(both), 1)))
    # no need to handle the case where only post-operative truth exists and the pre-operative data lies outside this batch
    out.append('')
    out.append('== per-patient, per-series mapping detail ==')
    rows = list(csv.DictReader(open(os.path.join(BASE, 'scripts', 'dicom_inventory.csv'), encoding='utf-8-sig')))
    byp = defaultdict(list)
    for r in rows:
        byp[r['patient']].append(r['desc'])
    for pat in sorted(byp):
        out.append('-- %s [%s]' % (pat, pre[pat]['vendor']))
        for d, c in Counter(byp[pat]).most_common():
            out.append('     %-42s x%-3d -> %s' % (d, c, vendor_of(d)))
    open(os.path.join(BASE, 'vendor_map2.txt'), 'w', encoding='utf-8').write('\n'.join(out))
    json.dump({k: v['vendor'] for k, v in pre.items()},
              open(os.path.join(BASE, 'vendor_map.json'), 'w', encoding='utf-8'),
              ensure_ascii=False, indent=1)
    print('agree', agree, 'disagree', dis)

if __name__ == '__main__':
    main()
