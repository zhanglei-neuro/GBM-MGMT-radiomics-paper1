# -*- coding: utf-8 -*-
r"""
patient_ids.py — 患者标识归一化

公开数据（Zenodo 数据集与期刊补充表）使用去标识编号 **P001–P109**
（见 `zenodo/README_DATASET.md` 与补充表 `TABLE_LEGENDS.md`）。
本地工作环境的数据表使用姓名拼音。把「编号 ↔ 本地姓名」的对照放在

    github/tools/patient_ids_local.csv      列: patient_id_public,patient_name_original

该文件含 PHI，**不要提交**（已列入 `github/.gitignore`）。
若本地数据本身已是编号（例如直接使用 Zenodo 发布版），无需该文件。

用法
----
    from patient_ids import to_id, to_name, normalize, drop_uih

    to_id('P001')        -> 'P001'      已是编号，幂等
    to_id('<PATIENT>')   -> 'P001'      本地姓名 → 编号（不区分大小写）
    to_name('P001')      -> '<PATIENT>' 编号 → 本地姓名（构造本地路径用）

三条约定
--------
1. **数据键**（队列成员、逐例映射的键、图表选例）统一用 Pxxx：
   读入患者索引的数据表后先 `normalize()`，再参与后续逻辑。
   这样同一份代码对本地（姓名）与发布（编号）两种数据源都成立。
2. **文件路径**必须用本地目录名：DICOM / NIfTI 目录以姓名命名，
   所以构造路径时用 `to_name(pid)`。
3. **队列排除一律按 vendor 判断**（见 `drop_uih` / `UIH_MASK`），不得按姓名：
   发布版数据中姓名已替换为编号，按姓名排除会**静默失效**并改变 n。
"""
import os
import csv

_HERE = os.path.dirname(os.path.abspath(__file__))

UIH = 'UIH'          # United Imaging 的 vendor 标签（cohort_MGMT.csv 的 vendor 列取值）

_LOCAL_CANDIDATES = [
    os.environ.get('PAPER1_PATIENT_MAP', ''),
    os.path.join(_HERE, 'patient_ids_local.csv'),
    os.path.join(_HERE, '..', '..', '_internal', 'patient_ids_local.csv'),
]

_ALIAS_CANDIDATES = [
    os.environ.get('PAPER1_PATIENT_ALIAS', ''),
    os.path.join(_HERE, 'patient_alias_local.csv'),
    os.path.join(_HERE, '..', '..', '_internal', 'patient_alias_local.csv'),
]


def load_map():
    """返回 (id2name, name2id)；找不到对照文件时返回两个空字典。"""
    id2name, name2id = {}, {}
    for p in _LOCAL_CANDIDATES:
        if p and os.path.isfile(p):
            with open(p, encoding='utf-8-sig', newline='') as f:
                r = csv.reader(f)
                next(r, None)
                for row in r:
                    if len(row) >= 2 and row[0].strip() and row[1].strip():
                        pid, nm = row[0].strip(), row[1].strip()
                        id2name[pid.upper()] = nm
                        name2id[nm.lower()] = pid
            break
    return id2name, name2id


ID2NAME, NAME2ID = load_map()
HAS_MAP = bool(ID2NAME)


def to_id(x):
    """姓名 → Pxxx；已是 Pxxx 或未知标识则原样返回。"""
    s = str(x).strip()
    if not s:
        return x
    if s.upper() in ID2NAME:
        return s.upper()
    return NAME2ID.get(s.lower(), x)


def to_name(x):
    """Pxxx → 本地姓名；非编号或未知则原样返回。"""
    s = str(x).strip()
    return ID2NAME.get(s.upper(), x)


def load_alias():
    r"""
    临床表拼写变体 → Pxxx。

    个别病例在临床表 `GBM_clinical_data.xlsx` 的 Name 列中与影像队列目录名的拼写
    不一致（例如首个音节不同、少一个字母）。把这类对照放在可选文件

        github/tools/patient_alias_local.csv      列: clinical_name,patient_id_public

    该文件**不含 PHI**（只有拼写变体与已有编号），可随仓库发布；缺省为空。
    """
    alias = {}
    for p in _ALIAS_CANDIDATES:
        if p and os.path.isfile(p):
            with open(p, encoding='utf-8-sig', newline='') as f:
                r = csv.reader(f)
                next(r, None)
                for row in r:
                    if len(row) >= 2 and row[0].strip() and row[1].strip():
                        alias[row[0].strip()] = row[1].strip()
                        alias[row[0].strip().lower()] = row[1].strip()
            break
    return alias


def resolve(x):
    """先查拼写变体表，再走 to_id；返回 Pxxx（未知标识原样返回）。"""
    s = str(x).strip()
    if not s:
        return x
    a = load_alias()
    return to_id(a.get(s, a.get(s.lower(), s)))


def normalize(df, col=None, inplace=True):
    """把 DataFrame 索引（或 `col` 列）中的患者标识统一为 Pxxx。"""
    out = df if inplace else df.copy()
    if col is None:
        out.index = [to_id(p) for p in out.index]
    else:
        out[col] = [to_id(p) for p in out[col]]
    return out


def norm_seq(seq):
    """把任意患者标识序列统一为 Pxxx。"""
    return [to_id(x) for x in seq]


def drop_uih(df, col='vendor'):
    """剔除单例 United Imaging 病例（按厂商列判断，不按姓名）。"""
    if df is None or col not in getattr(df, 'columns', []):
        return df
    return df[df[col].astype(str).str.strip() != UIH].copy()


def uih_mask(ids, co, col='vendor'):
    """给定患者标识序列与以 patient 为索引的队列表，返回 UIH 的布尔列表。"""
    return [str(co.loc[p, col]).strip() == UIH for p in ids]
