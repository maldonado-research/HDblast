"""Short, floating-point PDE controls for compensated initial profiles.

No matter extension is evolved. Frozen continuum equations/finite differences
are reused; every background coefficient is reset from the same zero-bump
transport formulation used to create the seed. No old payload is changed.
"""
from pathlib import Path
import argparse, hashlib, json, math, sys, time
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'frozen'))
import registered_solver as R
import balanced_constraint_seed as B

def profile(epsilon, rtol=8e-14):
    sol,meta=B.run(epsilon,rtol=rtol,transport=True)
    zb=sol.y[4,-1]
    def sample(z):
        targets=np.asarray(z)+zb
        if np.min(targets)<sol.y[4,0]:raise ValueError('Grid extends beyond cone start')
        indices=np.clip(np.searchsorted(sol.y[4],targets)-1,0,len(sol.t)-2)
        ys=np.array([brentq(lambda y:sol.sol(y)[4]-target,sol.t[i],sol.t[i+1],xtol=2e-14)
                     for target,i in zip(targets,indices)])
        return sol.sol(ys)
    return sample,meta

def setup(epsilon,hmin,L,stretch=.05):
    s=R.Solver(.001,hmin,L,stretch,dc=0.,ko=0.)
    reference,zero_meta=profile(0.)
    r,ry,psi,sy,zz,I=reference(s.z)
    s.bg=dict(rho=r,phi=psi-1,psi=psi,phiz=r*sy,Hc=ry,lnrho=np.log(r))
    s.rho=r;s.hc=ry;s.phz=r*sy;s.phi=psi-1;s.rb=r[-1];s.pb=s.phi[-1];s.rho2=r*r
    s.pot=R.derivatives(s.phi);s.s0=2*R.W(s.pb)+.001*(1+R.C*s.pb)
    s.s10=2*(s.pb*s.pb-1)+.001*R.C;s.s20=4*s.pb
    s.meta=zero_meta
    if epsilon:
        sample,seed_meta=profile(epsilon);rs,vs,ps,ss,zs,Is=sample(s.z)
    else:
        rs,vs,ps,ss,zs,Is=r.copy(),ry.copy(),psi.copy(),sy.copy(),zz.copy(),I.copy();seed_meta=zero_meta
    raw=np.zeros((6,s.n));raw[0]=np.log1p((rs-r)/r);raw[1]=raw[0];raw[2]=ps-psi
    state=raw.copy()
    x=np.clip((s.z+.99*L)/(.14*L),0,1);taper=x**3*(10-15*x+6*x*x)
    state[:3]*=taper;state[:,:2]=0
    ga,gf,_,_=s.bc(state[1,-1],state[2,-1],0.,0.)
    geometric_slopes=np.array([vs[-1]-ry[-1],rs[-1]*ss[-1]-r[-1]*sy[-1]])
    meta=dict(reference=zero_meta,seed=seed_meta,
              initialization='Dense ODE interpolation; separate inversion at identical shell-anchored conformal z; log1p warp and shifted scalar differences.',
              analytic_boundary_slope_minus_target=(geometric_slopes-[ga,gf]).tolist(),
              initial_shell_scalar_displacement=float(raw[2,-1]),
              initial_max_scalar_displacement=float(np.max(np.abs(raw[2]))),
              taper='C2 quintic smoothstep on [-.99L,-.85L]; modifies data outside the shell causal window.',
              earliest_taper_influence_on_shell=.85*L,
              initial_acceleration_bulk_max=float(np.max(np.abs(s.rhs(state)[3:,s.z>-.2]))))
    return s,state,meta

def diagnostics(s,v,t):
    row,H,M,C=s.diagnostics(v,t)
    scale=1+6*s.hc*s.hc+s.phz*s.phz
    core=s.z>-.2
    row.update(core_H=float(np.max(np.abs(H[core])/scale[core])),
               core_M=float(np.max(np.abs(M[core])/scale[core])),
               core_raw_H=float(np.max(np.abs(H[core]))),core_raw_M=float(np.max(np.abs(M[core]))),
               max_scalar_displacement=float(np.max(np.abs(v[2]))))
    a,b,f,pa,pb,pf=v
    ga,gf,_,_=s.bc(b[-1],f[-1],pb[-1],pf[-1])
    phiz=s.phz+s.D1@f+s.g1*gf
    peak=int(np.argmax(abs(phiz)))
    row['peak_scalar_gradient_z']=float(s.z[peak]);row['peak_scalar_gradient']=float(phiz[peak])
    row['shell_phi_velocity_over_H0']=float(pf[-1]*np.exp(-b[-1]))
    return row,H,M,C

def main():
    p=argparse.ArgumentParser();p.add_argument('--epsilon',type=float,default=1e-4)
    p.add_argument('--hmin',type=float,default=2e-4);p.add_argument('--L',type=float,default=6.)
    p.add_argument('--tf',type=float,default=2.5);p.add_argument('--output',required=True)
    p.add_argument('--init-only',action='store_true');p.add_argument('--max-f',type=float,default=.1)
    args=p.parse_args();out=Path(args.output);out.parent.mkdir(parents=True,exist_ok=True)
    s,v,meta=setup(args.epsilon,args.hmin,args.L)
    if args.tf>=.85*args.L:raise ValueError('Requested time reaches taper causal window')
    initial=v.copy();diag,H,M,C=diagnostics(s,v,0.)
    print(json.dumps({'initial':diag,'seed':meta}),flush=True)
    dt0=.4*s.zp.min();steps=0 if args.init_only else math.ceil(args.tf/dt0)
    dt=args.tf/steps if steps else 0.;every=max(1,round(.025/dt)) if dt else 1
    rows=[];snapshots={};tau=0.;start=time.monotonic();reason='init_only' if args.init_only else 'tf'
    for i in range(steps+1):
        t=i*dt
        if i%every==0 or i==steps:
            diag,H,M,C=diagnostics(s,v,t)
            row=dict(time=t,H0tau=tau,phi_b=s.pb+v[2,-1],scalar_deviation=v[2,-1],
                     H_over_H0=(1+v[3,-1])*np.exp(-v[1,-1]),**diag);rows.append(row)
            if len(rows)%20==1:print(json.dumps({'progress':row,'seconds':time.monotonic()-start}),flush=True)
            if len(rows)%20==1 or i==steps:
                snapshots['state_'+str(i)]=v.copy();snapshots['constraints_'+str(i)]=np.array([H,M,*C])
            if not np.isfinite(v).all() or np.max(np.abs(v[2]))>args.max_f:
                reason='amplitude_or_finite_guard';break
        if i==steps:break
        tau0=np.exp(v[1,-1]);k1=s.rhs(v);k2=s.rhs(v+.5*dt*k1);k3=s.rhs(v+.5*dt*k2);k4=s.rhs(v+dt*k3)
        v+=dt*(k1+2*k2+2*k3+k4)/6;tau+=dt*.5*(tau0+np.exp(v[1,-1]))
    result=dict(status='FLOATING_COMPENSATED_SEED_EVOLUTION',epsilon=args.epsilon,hmin=args.hmin,L=args.L,
                dt=dt,nodes=s.n,metadata=meta,stop_reason=reason,records=rows,runtime_seconds=time.monotonic()-start,
                scope='Classical matter-free short-time numerical control; no interval certificate, full nonlinear fate, or reheating.',
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(str(out)+'.npz',z=s.z,background=np.array([s.rho,s.phi,s.hc,s.phz]),initial=initial,final=v,**snapshots)
    print(json.dumps({'stop_reason':reason,'last':rows[-1],'runtime_seconds':result['runtime_seconds']}),flush=True)

if __name__=='__main__':main()
