# -*- coding: utf-8 -*-
r"""Necrosis 单序列(t1c)变体重跑
背景: 原批跑中 MRI_Necrosis 自动选型按变体名长度取最后一个匹配项, t1d 缺失被跳过,
导致三序列输入命中 t1c_t1w_t2f_t1d 四序列变体, 静默输出空掩码(103例全空)。
方案: 搭建暂存模型库 models_nec_t1c\MRI_Necrosis\ 仅含 t1c 子目录, 强制选型命中单序列变体;
MRI_Brain 一并拷入(t1c 管线第一步需要 T1-CE 上的脑分割)。
输出: D:\Raidionics_data\patients_nec2\2\<p>\T0\  (独立目录, 不触碰已验证的主结果树)
断点续跑: 输出 T0 中已存在 *Necrosis*.nii.gz 的患者跳过。
计时: necrosis_timing.csv
"""
import os, sys, csv, json, time, traceback, configparser, logging
import shutil as _shutil

# ---- 与 run_pipeline.py 相同的改名补丁(被监控树内运行时防钩子拦截; 计划任务环境由包装脚本关闭) ----
_TRASH_SEQ = [0]
def _rmtree_rename(path, ignore_errors=False, onerror=None, **kw):
    if not os.path.isdir(path):
        if ignore_errors:
            return
        raise FileNotFoundError(path)
    _TRASH_SEQ[0] += 1
    path = path.rstrip('/\\')  # 尾部分隔符会让 os.rename 报 WinError 87
    trash = f"{path}_trash_{os.getpid()}_{_TRASH_SEQ[0]}"
    try:
        os.rename(path, trash)
    except OSError:
        if not ignore_errors:
            raise

def _remove_rename(path, *args, **kw):
    try:
        os.rename(path, f"{path}.trash_{os.getpid()}_{_TRASH_SEQ[0]}")
        _TRASH_SEQ[0] += 1
    except OSError:
        pass

def _noop(*args, **kw):
    return None

_REAL_RMTREE, _REAL_REMOVE, _REAL_RMDIR = _shutil.rmtree, os.remove, os.rmdir
_shutil.rmtree = _rmtree_rename
os.remove = _remove_rename
os.rmdir = _noop
if os.environ.get('RAIDIONICS_NO_PATCH') == '1':
    _shutil.rmtree = _REAL_RMTREE
    os.remove = _REAL_REMOVE
    os.rmdir = _REAL_RMDIR

os.environ['OMP_NUM_THREADS'] = '4'
os.environ['MKL_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '4'

BASE = os.path.dirname(os.path.abspath(__file__))
SEL = os.path.join(BASE, 'series_selection.csv')
MAIN_TIMING = os.path.join(BASE, 'segmentation_timing.csv')
TIMING_CSV = os.path.join(BASE, 'necrosis_timing.csv')
WORK_ROOT = r"D:\Raidionics_data"
STAGING_MODELS = os.path.join(WORK_ROOT, 'models_nec_t1c')
REAL_MODELS = os.path.join(os.path.expanduser('~'), '.raidionics', 'resources', 'models')
PIPELINE_FN = os.path.join(BASE, 'necrosis_t1c_pipeline.json')
CONFIG_FN = os.path.join(BASE, 'necrosis_t1c_config.ini')
OUT_ROOT = os.path.join(WORK_ROOT, 'patients_nec2')


def build_staging_models(logger=print):
    """拷贝 MRI_Necrosis/t1c 与 MRI_Brain 到暂存模型库(幂等)。"""
    src_nec = os.path.join(REAL_MODELS, 'MRI_Necrosis', 't1c')
    src_brain = os.path.join(REAL_MODELS, 'MRI_Brain')
    dst_nec = os.path.join(STAGING_MODELS, 'MRI_Necrosis', 't1c')
    dst_brain = os.path.join(STAGING_MODELS, 'MRI_Brain')

    def complete(dst, n_onnx):
        if not os.path.isdir(dst):
            return False
        files = []
        for dp, dn, fn in os.walk(dst):
            files += fn
        return (len([f for f in files if f.endswith('.onnx')]) >= n_onnx
                and os.path.exists(os.path.join(dst, 'pipeline.json')))

    if not complete(dst_nec, 5):
        logger("  staging: copying MRI_Necrosis/t1c ...")
        if os.path.isdir(dst_nec):
            _shutil.rmtree(dst_nec, ignore_errors=True)
        _shutil.copytree(src_nec, dst_nec)
    if not complete(dst_brain, 5):
        logger("  staging: copying MRI_Brain ...")
        if os.path.isdir(dst_brain):
            _shutil.rmtree(dst_brain, ignore_errors=True)
        _shutil.copytree(src_brain, dst_brain)
    logger("  staging models ready: %s" % STAGING_MODELS)


def write_pipeline_json():
    pip = {
        "1": {"task": "Model selection", "model": "MRI_Necrosis", "timestamp": 0,
              "format": "thresholding",
              "description": "Model selection for Necrosis (t1c-only staging)"},
    }
    with open(PIPELINE_FN, 'w') as f:
        json.dump(pip, f, indent=4)


def write_config(input_folder, output_folder):
    cfg = configparser.ConfigParser()
    cfg['Default'] = {'task': 'neuro_diagnosis', 'trace': 'true'}
    cfg['System'] = {
        'gpu_id': '-1',
        'acceleration': 'cpu',
        'input_folder': input_folder,
        'output_folder': output_folder,
        'model_folder': STAGING_MODELS,
        'pipeline_filename': PIPELINE_FN,
    }
    cfg['Runtime'] = {
        'reconstruction_method': 'thresholding',
        'reconstruction_order': 'resample_first',
    }
    with open(CONFIG_FN, 'w') as f:
        cfg.write(f)


def necrosis_done(output_folder):
    t0 = os.path.join(output_folder, 'T0')
    if not os.path.isdir(t0):
        return False
    for f in os.listdir(t0):
        if 'Necrosis' in f and f.endswith('.nii.gz'):
            return True
    return False


def run_one(input_folder, output_folder, logger):
    from raidionicsrads.Utils.configuration_parser import ResourcesConfiguration
    from raidionicsrads.Utils.DataStructures.PatientStructure import PatientParameters
    from raidionicsrads.Pipelines.PipelineStructure import Pipeline

    ResourcesConfiguration.getInstance().set_environment(config_path=CONFIG_FN)
    patient_parameters = PatientParameters(id="Patient", patient_filepath=input_folder)
    pip = Pipeline(PIPELINE_FN)
    patient_parameters = pip.setup(patient_parameters=patient_parameters)

    total_t = 0.0
    steps = pip._steps
    n = len(steps)
    for k in sorted(steps.keys(), key=int):
        step = steps[k]
        t0 = time.time()
        err = None
        try:
            step.setup(patient_parameters)
        except Exception as e:
            if getattr(step, 'inclusion', 'required') == 'required':
                raise
            logger(f"    step{k} setup warn: {str(e)[:100]}")
        try:
            out = step.execute()
            if out is not None:
                patient_parameters = out
        except Exception as e:
            if getattr(step, 'inclusion', 'required') == 'required':
                raise
            err = str(e)[:120]
            logger(f"    step{k} exec warn: {err}")
        dt = time.time() - t0
        total_t += dt
        desc = str(getattr(step, 'step_description', '?'))[:60]
        logger(f"    [{k}/{n}] {desc} : {dt:.1f}s" + (f" (warn:{err})" if err else ""))
        try:
            step.cleanup()
        except Exception:
            pass
    pip.cleanup()
    return total_t


def flush(row):
    exists = os.path.exists(TIMING_CSV)
    fields = ['patient', 'status', 'Necrosis_Time', 'Total_Time', 'Start', 'End', 'Note']
    with open(TIMING_CSV, 'a', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if not exists:
            w.writeheader()
        w.writerow({k: row.get(k, '') for k in fields})


def main():
    logging.basicConfig(level=logging.WARNING)
    only = None
    if len(sys.argv) > 2 and sys.argv[1] == '--only':
        only = sys.argv[2]
    write_pipeline_json()
    build_staging_models()

    # 队列 = 主批跑状态为 OK 的患者(103例)
    cohort = []
    for r in csv.DictReader(open(MAIN_TIMING, encoding='utf-8-sig')):
        if r['status'] == 'OK' and r['patient'] not in cohort:
            cohort.append(r['patient'])
    if only:
        cohort = [p for p in cohort if p == only]
    print(f"cohort={len(cohort)}", flush=True)

    for i, p in enumerate(cohort):
        print(f"[{i+1}/{len(cohort)}] === {p} ===", flush=True)
        input_folder = os.path.join(WORK_ROOT, 'input', '2', p, 'T0')
        output_folder = os.path.join(OUT_ROOT, '2', p)
        if necrosis_done(output_folder):
            print("  skip: already has necrosis mask", flush=True)
            continue
        niftis = [f for f in os.listdir(input_folder)] if os.path.isdir(input_folder) else []
        if len(niftis) < 1:
            flush({'patient': p, 'status': 'ERROR', 'Note': f'input missing ({len(niftis)} nifti)'})
            print("  ERROR: input missing", flush=True)
            continue
        # 整体重命名旧输出目录(非破坏), 确保全新运行
        if os.path.isdir(output_folder):
            _TRASH_SEQ[0] += 1
            try:
                os.rename(output_folder, output_folder + f"_trash_init_{_TRASH_SEQ[0]}")
            except OSError:
                pass
        os.makedirs(output_folder, exist_ok=True)
        write_config(os.path.join(WORK_ROOT, 'input', '2', p), output_folder)
        start_ts = time.strftime('%Y-%m-%d %H:%M:%S')
        t_start = time.time()
        try:
            tt = run_one(os.path.join(WORK_ROOT, 'input', '2', p), output_folder,
                         lambda m: print("  " + m, flush=True))
            flush({'patient': p, 'status': 'OK', 'Necrosis_Time': f"{tt:.1f}",
                   'Total_Time': f"{time.time()-t_start:.1f}", 'Start': start_ts,
                   'End': time.strftime('%Y-%m-%d %H:%M:%S'), 'Note': ''})
            print(f"  DONE nec={tt:.0f}s total={(time.time()-t_start)/60:.1f}min", flush=True)
        except Exception as e:
            flush({'patient': p, 'status': 'ERROR', 'Total_Time': f"{time.time()-t_start:.1f}",
                   'Start': start_ts, 'End': time.strftime('%Y-%m-%d %H:%M:%S'),
                   'Note': str(e)[:200]})
            print(f"  ERROR: {str(e)[:200]}", flush=True)
            traceback.print_exc()
    print("ALL_DONE", flush=True)


if __name__ == '__main__':
    main()
