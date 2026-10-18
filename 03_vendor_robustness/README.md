# 厂商鲁棒性 (Vendor robustness)

GE 亚组 (n=32) 的逐折诊断：训练/测试类别分布、逐折 AUC 与方向、内层筛选方向的稳定性；以及 100 次标签置换零分布，用于判断亚组内低于随机水平的 AUC 是否可由随机波动解释。

## 本目录文件

- `ge_subgroup_fold_diagnostics.py`
- `ge_subgroup_permutation_null.py`
