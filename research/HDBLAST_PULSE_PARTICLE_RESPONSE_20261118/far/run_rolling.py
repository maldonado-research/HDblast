"""Paired original B1 matter runs forked on one common, pre-crossing numerical slice."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
B1='research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/'
PINS={B1+'evolve_b1.py':'a0ae79e12ff90f549a0d1ca674fa3285b8ad456b780880df405b9597d3434f0f',
      B1+'static_w.py':'ac4f5df5d8a1fcedfbd9503caf4717c1473fde6699388e095618f269b6577b1f',
      B1+'shell_matter.py':'58c4c4cc18132f86aaa45ac9f0930490ef3debe160fe8f6ed59274451aadc7d1',
      'research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/rolloff5d_v1.py':'9ea4d8d6ef36ce3d242fbb0e09f27a944b571fdfb8646e3d71195b60b089080b'}
# The producer's RK4 equations, matter-event handler and three junctions are unmodified.
MARK='        n = 0; last_log = time.time()'
RESTORE="""
        if self.branch != 'capture':
            with np.load(self.checkpoint, allow_pickle=False) as cp:
                Yv = cp['Yv'].copy(); m=cp['m'].copy(); T=float(cp['T']); tau=float(cp['tau'])
                n=int(cp['step']); nsnap=int(cp['nsnap']); A0=float(cp['A0'])
                hist.extend(cp['hist'].tolist()); self.mode=str(cp['mode'].item())
            if self.branch!='quiet':
                from rolling_data import add_pulse
                Yv, receipt=add_pulse(self,T,Yv,m,float(self.branch),self.response_options)
                Path(out_dir).mkdir(parents=True,exist_ok=True)
                (Path(out_dir)/'initial.json').write_text(json.dumps(receipt,indent=2)+'\\n')
"""
STEP='        while True:\n            if n % every == 0:'
CAPTURE="""
            if self.branch=='capture' and n==round(1.11/dt):
                if len(getattr(self.matter,'events',())) or m.size!=1:
                    raise ValueError('Expected preproduction state')
                np.savez(self.checkpoint,Yv=Yv,m=m,T=T,tau=tau,step=n,nsnap=nsnap,A0=A0,
                         hist=np.array(hist),mode=self.mode)
                return
"""
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('repo',type=Path);p.add_argument('workspace',type=Path)
    p.add_argument('--grid',choices=['coarse','fine'],required=True)
    p.add_argument('--branch',choices=['capture','quiet','0.001','-0.001','0.0005','-0.0005'],required=True)
    p.add_argument('--end',type=float,default=1.62)
    p.add_argument('--kind',choices=['incident','velocity','metric','metric_field'],default='incident')
    p.add_argument('--control',type=float,default=-.85)
    p.add_argument('--spread',type=float,default=.25)
    p.add_argument('--label')
    args=p.parse_args()
    work=args.workspace.resolve()/args.grid
    work.mkdir(parents=True,exist_ok=True)
    program=work/'solver'
    if args.branch=='capture':
        program.mkdir(exist_ok=True)
        for name,digest in PINS.items():
            data=(args.repo/name).read_bytes()
            if hashlib.sha256(data).hexdigest()!=digest:raise ValueError(name)
            if name.endswith('/evolve_b1.py'):
                text=data.decode()
                if text.count(MARK)!=1 or text.count(STEP)!=1:raise ValueError('Unknown B1 source')
                text=text.replace(MARK,MARK+RESTORE).replace(STEP,'        while True:'+CAPTURE+'            if n % every == 0:')
                data=text.encode()
            (program/Path(name).name).write_bytes(data)
    else:
        # Inject the numerical construction only in the independent branch process.
        if not (work/'checkpoint.npz').exists():raise FileNotFoundError('Run capture first')
    driver="""import sys,math,json
from pathlib import Path
from evolve_b1 import Run
from shell_matter import ChiGas
sys.path.insert(0,sys.argv[1])
base=Path(sys.argv[2]); name=sys.argv[3]; fine=(sys.argv[4]=='fine')
run=Run(delta=.1,c=.5975949350280132,d=-3.106933495673783,dc=.01,
        dzf=.0005 if fine else .001,dzc=.002 if fine else .004,L=16.,
        kappa=10.,project=15.9 if name=='capture' else None,cfl=.5,order=4,xc=1.9,
        matter=lambda rb:ChiGas(100.,.5,1.,6.49e-6,rb),log=lambda msg: None)
run.branch=name;run.checkpoint=base/'checkpoint.npz';run.response_options={'kind':sys.argv[6],'control':float(sys.argv[7]),'spread':float(sys.argv[8])}
run.evolve(T_final=float(sys.argv[5]),rec_dT=.001,tag='rolling',out_dir=str(base/sys.argv[9]),wall_limit=1000)
"""
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')
    subprocess.run([sys.executable,'-B','-c',driver,str(HERE),str(work),args.branch,args.grid,str(args.end),args.kind,str(args.control),str(args.spread),args.label or args.branch],cwd=program,check=True,env=env)

if __name__=='__main__':main()
