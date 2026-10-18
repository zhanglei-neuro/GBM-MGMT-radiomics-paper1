# 工具 (Tools)

`patient_ids.py` 患者标识归一化（编号 ↔ 本地姓名、按厂商排除单例病例）；
`check_phi.py` 上传前患者隐私自检；`sanitize_phi.py` 姓名占位符替换；
`configure_paths.py` 批量改写硬编码路径。

**每次 `git push` 之前请先运行 `check_phi.py`。**

## 本目录文件

- `patient_ids.py` —— 被 01/02/03/06/07/08/09 共 18 个脚本导入，见仓库 README §4
- `patient_alias_local.csv.example` —— 临床表拼写变体表模板（可选）
- `check_phi.py`
- `configure_paths.py`
- `sanitize_phi.py`

## 本地运行需要的对照表（**不要提交**）

| 文件 | 内容 | 获取方式 |
|---|---|---|
| `patient_ids_local.csv` | 两列 `patient_id_public,patient_name_original` | 从你的安全位置复制到此，或设 `PAPER1_PATIENT_MAP` 指向它 |
| `patient_alias_local.csv` | 两列 `clinical_name,patient_id_public`（临床表与队列目录名拼写不一致者） | 同上，或设 `PAPER1_PATIENT_ALIAS` |

两者都已列入 `.gitignore`。缺少第一个文件时脚本会退化为"标识原样使用"，
即只能读取已使用 `P001`–`P109` 编号的发布版数据。
