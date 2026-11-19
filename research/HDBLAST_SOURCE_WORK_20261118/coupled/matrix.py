"""Local run scheduler for the published set of comparisons."""
import argparse
import concurrent.futures
from pathlib import Path
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repo',type=Path)
    parser.add_argument('work',type=Path)
    parser.add_argument('--jobs',type=int,default=4)
    parser.add_argument('--sensitivity',action='store_true')
    opt=parser.parse_args()
    work=opt.work.resolve()
    driver=Path(__file__).with_name('run_coupled.py')
    def initial(grid):
        dest=work/grid/'initial_T0.700.npz'
        if not dest.exists():
            subprocess.run([sys.executable,'-B',driver,opt.repo,work,'--grid',grid,'--branch','initial'],check=True)
    def evolve(row):
        grid,y,on,nk,branch=row
        folder=work/grid/f'k{nk}_y{y}_b{on}_t0.700';folder.mkdir(exist_ok=True)
        cmd=[sys.executable,'-B',driver,opt.repo,work,'--grid',grid,'--decay',str(y),'--nodes',str(nk),'--branch',branch]
        if not on:cmd.append('--passive')
        completed= folder/('prepare.json' if branch=='fork' else branch+'/wall_summary.json')
        if completed.exists():return
        with (folder/(branch+'.log')).open('w') as f:
            subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        print(*row,'complete',flush=True)
    if opt.sensitivity:
        groups=[('fine',48,.70),('coarse',24,.65),('coarse',24,.75)]
        # Both coarse preparation times share the native solver's static cache.
        # Populate each grid on one process before running independent branches.
        for grid in ('coarse','fine'):
            initial(grid)
        def sensitivity(c):
            grid,nodes,prep=c
            folder=work/grid/f'k{nodes}_y1_b1_t{prep:.3f}'
            folder.mkdir(exist_ok=True)
            for stage in ('initial','fork','quiet','plus','minus','far'):
                if stage=='initial' and (work/grid/f'initial_T{prep:.3f}.npz').exists():
                    continue
                cmd=[sys.executable,'-B',driver,opt.repo,work,'--grid',grid,'--decay','1','--nodes',str(nodes),'--prep-T',str(prep),'--branch',stage]
                with (folder/(stage+'.log')).open('w') as log:
                    subprocess.run(cmd,stdout=log,stderr=log,check=True)
                print(*c,stage,'complete',flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=opt.jobs) as pool:
            list(pool.map(sensitivity,groups))
        return
    with concurrent.futures.ThreadPoolExecutor(max_workers=opt.jobs) as pool:
        for _ in pool.map(initial,('coarse','fine')): pass
        # Submit each run as soon as its own common matter checkpoint is available.
        cases=[(g,y,on,24) for g in ('coarse','fine') for y,on in [(1,1),(1,0),(0,1),(0,0)]]
        tasks={pool.submit(evolve,(*c,'fork')):c for c in cases}
        while tasks:
            done,_=concurrent.futures.wait(tasks,return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in done:
                c=tasks.pop(fut);fut.result()
                if c:
                    for arm in ('quiet','plus','minus','far'):
                        tasks[pool.submit(evolve,(*c,arm))]=None


if __name__=='__main__':
    main()
