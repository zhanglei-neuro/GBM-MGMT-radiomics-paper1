# ComBat 调和 (ComBat harmonization)

> **本目录不含稿件所报告 ComBat 结果的实现。** 真正的实现在主分析管线内：
> **`../02_nested_cross_validation/paper1_nested_cv.py` 第 [7] 节 `ComBat harmonisation (neuroCombat)`**。

## 稿件所报告的方法（权威口径）

- 实现：**neuroCombat 0.2.12**（经验贝叶斯收缩；Johnson 2007 / Fortin 2017）
- 调用方式：每外折在**训练折内**拟合，再经 `neuroCombatFromTraining` 变换测试折 → 无参数泄漏
- 批次变量：仅传入厂商标签，**不传入任何生物学协变量**（MGMT 标签全程不参与调和）
- 报告内容：调和后嵌套 CV（logistic / XGBoost）、调和后 LOVO、调和后厂商显著特征数（对比调和前），
  以及转导式替代方案 **VWS**（厂商内 Z 标准化，均值/方差仅在训练折估计）的 LOVO 对照
- 产物：`zenodo/03_intermediate_results/metrics_combat.json`

复现命令：

```bash
python 02_nested_cross_validation/paper1_nested_cv.py
```

## 本目录文件

- `combat_protected_covariate_sensitivity.py` — **稿件 §3.7 报告的保护协变量敏感性分析**。
  在同一嵌套 CV 框架内比较两臂：A 仅传厂商批次变量（主分析设定）；B 在训练折内额外以 MGMT 为
  discrete covariate（`categorical_cols=['MGMT']`）。均在训练折拟合、经 `neuroCombatFromTraining`
  变换测试折。结论：保护协变量**未**恢复判别力（logistic ΔAUC −0.017，p = 0.77；
  XGBoost ΔAUC −0.127，p = 0.09）。结果文件
  `zenodo/03_intermediate_results/combat_protected_covariate_sensitivity.json`。
- `legacy_simplified_combat_exploratory.py` — **早期探索性实现，未用于稿件**。
  实现的是无经验贝叶斯收缩的简化 location/scale 版本（Johnson 2007），可 fit/transform，
  下游复用与主分析相同的嵌套 CV 管线。保留仅用于溯源；其输出（`predictions_combat.csv`、
  旧版 `metrics_combat.json`）**未被稿件采用**，请勿用于复现稿件结果。

> 说明：本目录原文件名为 `combat_harmonization.py`，容易与稿件报告的方法混淆，
> 2026-10-12 修订时改名并加注废弃横幅。
