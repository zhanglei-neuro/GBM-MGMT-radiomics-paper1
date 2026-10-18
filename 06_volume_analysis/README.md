# 体积生物标志物 (Volume biomarkers)

7 项预设体积标志物：TumorCore / Necrosis / Edema / FLAIRChanges 体积、总体肿瘤体积、总体肿瘤/水肿比、强化核心/水肿比。体积由掩膜在**原生网格上逐体素计数**得到，避免重采样带来的体积漂移。统计：Mann–Whitney、单变量 AUC + bootstrap CI、族内 BH-FDR。

## 本目录文件

- `make_per_case_volumes.py`
- `volume_biomarker_analysis.py`
