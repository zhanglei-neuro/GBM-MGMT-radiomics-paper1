# 嵌套交叉验证主分析 (Nested cross-validation)

论文主分析。全局标签无关预处理 → 外层分层 5 折 × 5 次重复 → 内层单变量 |AUC−0.5| top-200 → |Spearman|>0.9 贪心去相关（≤40 特征）→ 3 折超参搜索。全部预处理仅在训练折拟合。
模型：L2 逻辑回归 / XGBoost / 随机森林；对照：临床基线（年龄+性别）与联合模型。
指标：重复均值±SD、合并 OOF AUC、DeLong 检验、类别无关 NRI、IDI、校准斜率/截距、Brier、bootstrap 95%CI、LOVO。

Entry point: `paper1_nested_cv.py`

## 本目录文件

- `paper1_nested_cv.py`
