# -*- coding: utf-8 -*-
r"""Batch launcher: Necrosis rerun on t1c, with stdout/stderr also written to run_necrosis_batch.log.

Differs from pointing the scheduled task straight at necrosis_rerun.py:
  - sets RAIDIONICS_NO_PATCH=1, which turns off the "rename instead of delete" patch
    inside necrosis_rerun.py (the task scheduler is non-interactive, which avoids a
    pile-up of _trash_ directories);
  - caps the BLAS thread count at 4 to avoid thread contention in multi-process runs;
  - tees all output to disk so every case can be audited individually.
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
