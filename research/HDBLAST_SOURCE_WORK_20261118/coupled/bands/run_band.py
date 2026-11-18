"""Regenerate one wider finite-band experiment using the released native driver."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def native_call():
    # The native worker imports quantum_shell from HERE, the only substitution.
    src=HERE.parent/'run_coupled.py'
    spec=importlib.util.spec_from_file_location('native_driver',src)
    driver=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    driver.HERE=HERE
    sys.argv[1:]=sys.argv[2:]
    driver.main()


def regenerate():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('repo',type=Path)
    p.add_argument('work',type=Path)
    p.add_argument('--edge',choices=(12,14),type=int,required=True)
    args=p.parse_args()
    def run_grid(grid):
        base=args.work.resolve()/f'K{args.edge}'
        for stage in ('initial','fork','quiet','plus','minus','far'):
            cmd=[sys.executable,'-I','-B',__file__,'--native',args.repo.resolve(),base,
                 '--grid',grid,'--nodes','32','--decay','1','--branch',stage]
            env=dict(os.environ,HDBLAST_COMOVING_BAND_EDGE=str(args.edge))
            subprocess.run(cmd,env=env,check=True)
    with ThreadPoolExecutor(max_workers=2) as queue:
        list(queue.map(run_grid,['coarse','fine']))


if __name__=='__main__':
    if sys.argv[1:2]==['--native']:
        native_call()
    else:
        regenerate()
