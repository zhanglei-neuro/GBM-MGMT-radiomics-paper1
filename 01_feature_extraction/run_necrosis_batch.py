# -*- coding: utf-8 -*-
r"""批处理启动包装：Necrosis t1c 重跑，stdout/stderr 同时写入 run_necrosis_batch.log。

与直接将 necrosis_rerun.py 作为计划任务目标的区别：
  - 设置 RAIDIONICS_NO_PATCH=1，关闭 necrosis_rerun.py 内的"改名代替删除"补丁
    （计划任务环境无交互，避免 _trash_ 目录堆积）；
  - 限制 BLAS 线程数（4），避免多进程环境下的线程争抢；
  - 用 Tee 把全部输出落盘，便于逐例审计。
"""
import os
import sys

os.environ['RAIDIONICS_NO_PATCH'] = '1'
os.environ['OMP_NUM_THREADS'] = '4'
os.environ['MKL_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'
os.environ['NUMEXPR_NUM_THREADS'] = '4'
os.environ['PYTHONIOENCODING'] = 'utf-8'
os.environ['PYTHONUTF8'] = '1'

LOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_necrosis_batch.log')


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            try:
                s.write(data)
                s.flush()
            except Exception:
                pass

    def flush(self):
        for s in self.streams:
            try:
                s.flush()
            except Exception:
                pass


logf = open(LOG, 'a', encoding='utf-8', buffering=1)
sys.stdout = Tee(sys.__stdout__, logf)
sys.stderr = Tee(sys.__stderr__, logf)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import necrosis_rerun  # noqa: E402

necrosis_rerun.main()
