# 特征提取与数据准备 (Feature extraction & data preparation)

从原始术前 DICOM 到 PyRadiomics 特征矩阵的全过程：DICOM 清单扫描 → 序列角色判定 (T1c 依据 `ContrastBolusAgent = Gd-DTPA`) → DICOM→NIfTI 转换 → 队列锁定与掩码齐全度校验 → 4148 维特征提取（Python 3.9 / pyradiomics 3.1.0，固定 bin 宽 25，不做强度归一化）。

Entry point: `extract_features_pyradiomics.py`

## 本目录文件

- `build_cohort.py`
- `build_vendor_map.py`
- `convert_dicom_to_nifti.py`
- `extract_features_pyradiomics.py`
- `scan_dicom_inventory.py`
- `select_series_roles.py`

### Necrosis 通道重跑（主批跑后的必要补救）

主批跑中 Raidionics 的 `MRI_Necrosis` 模型自动选型按变体名长度取最后一个匹配项；当输入缺少 t1d 序列时，选型落到 `t1c_t1w_t2f_t1d` 四序列变体，**静默输出空掩膜**（全部 90 例 Necrosis 掩膜体素数为 0，且进程返回码为 0、日志无报错）。因此主结果树的 Necrosis 通道不可用，须按下述流程重跑。

- `necrosis_rerun.py` — 单序列 (t1c) 变体重跑主脚本。搭建仅含 `t1c` 子目录的暂存模型库 `models_nec_t1c/`，强制选型命中单序列变体；输出写入**独立目录** `patients_nec2/`，不触碰已验证的主结果树；支持断点续跑（输出 `T0/` 下已有 `*Necrosis*.nii.gz` 的病例跳过）。
  用法：`python necrosis_rerun.py [--only <patient>]`
- `necrosis_t1c_pipeline.json` — 单步管线定义（仅 Model selection，`MRI_Necrosis`，thresholding 重建），由脚本运行时读取。
- `necrosis_config_template.ini` — Raidionics 运行时配置模板（`gpu_id=-1`、`acceleration=cpu`、`reconstruction_order=resample_first`）。占位符 `<WORK_ROOT>` / `<SCRIPT_DIR>` 由运行环境替换；脚本会按病例自动生成实际配置。
- `run_necrosis_batch.py` — 批处理启动包装：关闭重跑脚本内的改名补丁、限制 BLAS 线程数、把全部输出 tee 到 `run_necrosis_batch.log`。
- `necrosis_check.py` — 空掩码核查工具。统计各亚区（TumorCore / Necrosis / Edema / FLAIRChanges）的标签取值与非零体素数，用于判定某一通道是否整体失效。
  用法：`python necrosis_check.py <segmentation_tree> [report.txt]`

**结果（n=90）**：重跑后 84/90 例检出坏死区（中位 4.4 mL），6 例为真实无坏死（非技术失败），与主批跑 90/90 全空形成对照。

> ⚠️ **数据保护说明**：`necrosis_rerun.py` 在导入时对 `shutil.rmtree` / `os.remove` / `os.rmdir` 打补丁，将删除操作改为**重命名**（`*_trash_<pid>_<seq>`），以避免 Raidionics 内部清理逻辑误删已验证数据。设置环境变量 `RAIDIONICS_NO_PATCH=1` 可恢复原生行为（`run_necrosis_batch.py` 即采用此模式）。该补丁只改名、不删除，但会在输出树留下 `_trash_*` 目录，可用后自行清理。
