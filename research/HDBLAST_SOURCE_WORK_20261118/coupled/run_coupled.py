"""Run the original bulk RK4 with a Gaussian shell-matter sector and a common pulse fork."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import tempfile
import subprocess
import sys

HERE = Path(__file__).resolve().parent
OLD = Path("research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical")
SOURCE = Path("research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source")
PINS = {
    SOURCE / 'evolve_b1.py': 'a0ae79e12ff90f549a0d1ca674fa3285b8ad456b780880df405b9597d3434f0f',
    SOURCE / 'static_w.py': 'ac4f5df5d8a1fcedfbd9503caf4717c1473fde6699388e095618f269b6577b1f',
    SOURCE / 'shell_matter.py': '58c4c4cc18132f86aaa45ac9f0930490ef3debe160fe8f6ed59274451aadc7d1',
    Path("research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/rolloff5d_v1.py"):
        '9ea4d8d6ef36ce3d242fbb0e09f27a944b571fdfb8646e3d71195b60b089080b',
}
RESTORE = """
        if self.branch != 'initial':
            with np.load(self.fork,allow_pickle=False) as cp:
                Yv=cp['bulk']; m=cp['matter']; T=float(cp['T']); tau=float(cp['tau'])
                n=int(cp['step']); nsnap=int(cp['nsnap']); A0=float(cp['A0']); self.mode=str(cp['mode'].item())
            if self.branch=='fork':
                Ab,Bb,phib,AT,vT = self.shell_vals(T,Yv)
                m=self.matter.prepare(Ab,phib)
                rho,p,j=self.matter.totals(m,Ab,phib)
                dm=self.matter.rhs(m,Ab,phib,AT,vT,math.exp(Bb))
                dr=self.matter.totals_T(m,dm,Ab,phib,AT,vT,math.exp(Bb))
                self.preparation=dict(T=T,s=tau,energy_density=rho,pressure=p,force=j,derivatives=list(dr),
                                      phi=phib,H_over_H0=AT*self.rb/math.exp(Bb))
            elif self.branch!='quiet':
                from rolling_data import add_pulse
                amplitude=-.001 if self.branch=='minus' else .001
                options={'control':-.85,'spread':.25} if self.branch=='far' else {}
                Yv,self.match=add_pulse(self,T,Yv,m,amplitude,options)
"""
SAVE = """
            if self.branch in ('initial','fork') and n==round((self.prep_T if self.branch=='initial' else 1.11)/dt):
                np.savez(self.save_fork,bulk=Yv,matter=m,T=T,tau=tau,step=n,nsnap=nsnap,A0=A0,mode=self.mode)
                if self.branch=='fork':
                    np.savez_compressed(str(self.save_fork)+'.history.npz', cols=list(rows[-1].keys()), rec=np.array([[r.get(k,np.nan) for k in rows[-1]] for r in rows],dtype=str))
                return
"""


def patched(text):
    changes=[
        ("        max_mchi = 0.0",RESTORE+"        max_mchi = 0.0"),
        ("        while True:\n            if n % every == 0:","        while True:"+SAVE+"            if n % every == 0:"),
        ("R = float(m[0]); X = rho - R","R = self.matter.mult*float(m[0]) if self.matter.feedback else 0.; X = rho - R"),
        ("R_T=float(dm[0])","R_T=self.matter.mult*float(dm[0]) if self.matter.feedback else 0."),
        ("for k_ in ('n_hat', 'm_hat', 'rho_hat', 'p_hat', 'j_hat', 'jprod_hat', 'Q_hat', 'P_hat', 'ncoh', 'u_last'):","for k_ in ex:"),
    ]
    for old,new in changes:
        if text.count(old)!=1:
            raise ValueError("Unexpected original bulk source: "+old)
        text=text.replace(old,new)
    return text


def install(path, content):
    if path.exists() and path.read_bytes()==content:
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,delete=False) as temporary:
        temporary.write(content)
    os.replace(temporary.name,path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repo',type=Path)
    parser.add_argument('work',type=Path)
    parser.add_argument('--grid',choices=('coarse','fine'),required=True)
    parser.add_argument('--branch',choices=('initial','fork','quiet','plus','minus','far'),required=True)
    parser.add_argument('--nodes',type=int,default=24)
    parser.add_argument('--decay',type=float,choices=(0.,1.),default=0.)
    parser.add_argument('--passive',action='store_true')
    parser.add_argument('--prep-T',type=float,default=.70)
    cfg=parser.parse_args()
    repo,work=cfg.repo.resolve(),cfg.work.resolve()
    solver=work/cfg.grid/'solver'
    solver.mkdir(exist_ok=True,parents=True)
    branch=cfg.branch
    data=work/cfg.grid/f'k{cfg.nodes}_y{int(cfg.decay)}_b{int(not cfg.passive)}_t{cfg.prep_T:.3f}'
    data.mkdir(exist_ok=True)
    initial=solver.parent/f'initial_T{cfg.prep_T:.3f}.npz'
    fork=data/'fork.npz'
    if branch=='initial':
        for name,pin in PINS.items():
            content=(repo/name).read_bytes()
            if hashlib.sha256(content).hexdigest()!=pin:raise ValueError(name)
            if name.name=='evolve_b1.py':content=patched(content.decode()).encode()
            install(solver/name.name,content)
        orig=(repo/OLD/'rolling_data.py').read_bytes()
        if hashlib.sha256(orig).hexdigest()!=MATCH_PIN:raise ValueError('Incorrect pulse source')
        install(solver/'rolling_data.py',orig)
    event=json.loads((repo/OLD/'runs'/cfg.grid/'quiet/rolling_summary.json').read_bytes())['matter']['events'][0]
    config=dict(nodes=cfg.nodes,feedback=not cfg.passive,decay=cfg.decay,loga_ref=event['A_star'])
    env=dict(os.environ, OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    worker="""
import json,sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
from evolve_b1 import Run
from quantum_shell import QuantumShell
grid,branch=sys.argv[2:4]; initial,fork,root=map(Path,sys.argv[4:7]); cfg=json.loads(sys.argv[7])
run=Run(delta=.1,c=.5975949350280132,d=-3.106933495673783,dc=.01,dzf=.0005 if grid=='fine' else .001,
        dzc=.002 if grid=='fine' else .004,L=16.,kappa=10.,project=15.9 if branch=='initial' else None,
        cfl=.5,order=4,xc=1.9,matter=lambda rb:QuantumShell(rb,**cfg),log=lambda msg:None)
run.branch=branch; run.prep_T=float(sys.argv[8])
run.fork=initial if branch=='fork' else fork; run.save_fork=initial if branch=='initial' else fork
out=root/branch
run.evolve(T_final=1.62,rec_dT=.001,tag='wall',out_dir=str(out),wall_limit=2000)
if branch=='fork':
    (root/'prepare.json').write_text(json.dumps(run.preparation,indent=2)+'\\n')
if branch not in ('initial','fork'):
    if hasattr(run,'match'):(out/'match.json').write_text(json.dumps(run.match,indent=2)+'\\n')
    d=json.loads((out/'wall_summary.json').read_bytes())
    if d['stop_reason']!='T_final':raise ValueError(d['stop_reason'])
"""
    subprocess.run([sys.executable,'-B','-c',worker,str(HERE),cfg.grid,branch,str(initial),str(fork),
                    str(data),json.dumps(config),str(cfg.prep_T)],cwd=solver,check=True,env=env)


MATCH_PIN='2f024e981c6809f60147d8da9a639662240ae7643e0bb6b65b7feb67c3b3447e'
if __name__=='__main__':
    main()
