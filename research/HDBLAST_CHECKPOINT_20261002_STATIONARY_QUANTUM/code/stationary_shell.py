#!/usr/bin/env python3
"""Registered finite-amplitude stationary shell continuation (mathematical benchmark)."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]
DESITTER = 'research/HDBLAST_CHECKPOINT_20261002_DESITTER'
INPUT_FILES = [
    DESITTER+'/code/desitter_sources.py',
    DESITTER+'/sensitivity/compute_sensitivity.py',
    DESITTER+'/sensitivity/CHECKS.json',
    DESITTER+'/sensitivity/STATIC_SENSITIVITY.md',
    'research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/static_bridge/STATIC_DESITTER_SOURCE_BRIDGE.md',
    'research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/ACTION_AND_JUNCTION_SPECIFICATION.md',
]
LOCAL_FILES = ['REGISTRATION.md', 'EXPERIMENT.json', 'FIXED_MODEL.json', 'code/stationary_shell.py',
               'code/validate_stationary.py', 'requirements.txt']
EH0 = -1.3366923651588084e-23
Y0 = 30.276680395885563
ETA0 = -8.405268308670538e-5
Z0 = 5.924014794328956e-5
REFERENCE = '0.00011848029588657912'
DELTA = .001
DETUNING_C = .5975949350280132
SCALES = np.array([Z0, DELTA*DETUNING_C/2])

def utc():
    return datetime.now(timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False)+'\n')

def module(path, name):
    spec=importlib.util.spec_from_file_location(name, path)
    loaded=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded

def freeze(repo):
    target=ROOT/'PREREGISTRATION.json'
    if target.exists() or (ROOT/'outputs').exists():
        raise RuntimeError('Registration/output exists; refusing retrospective registration')
    doc={'registered_at_utc':utc(),
         'source_commit':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
         'local_sha256':{f:sha(ROOT/f) for f in LOCAL_FILES},
         'source_sha256':{f:sha(repo/f) for f in INPUT_FILES},
         'experiment':json.loads((ROOT/'EXPERIMENT.json').read_text()),
         'execution_rule':'Publish this registration and exact source before any new source evaluation or radial integration.'}
    write_json(target,doc)
    return doc

def verify_inputs(repo):
    doc=json.loads((ROOT/'PREREGISTRATION.json').read_text())
    for name,digest in doc['local_sha256'].items():
        if sha(ROOT/name)!=digest:
            raise RuntimeError('Frozen local input mismatch: '+name)
    for name,digest in doc['source_sha256'].items():
        if sha(repo/name)!=digest:
            raise RuntimeError('Inherited source input mismatch: '+name)
    return doc

class Model:
    def __init__(self, repo, dps=60, terms=128):
        self.bulk=module(repo/DESITTER/'sensitivity/compute_sensitivity.py','inherited_bulk')
        self.action=module(repo/DESITTER/'code/desitter_sources.py','inherited_action')
        self.archive=json.loads((repo/DESITTER/'sensitivity/CHECKS.json').read_text())['reference']
        self.dps=dps
        self.terms=terms
        self.source_cache={}

    def source(self, eta, z, b):
        key=(float(eta),float(z),int(b))
        if key in self.source_cache:
            return self.source_cache[key]
        with mp.workdps(self.dps):
            r=mp.mpf(REFERENCE)
            zz=mp.mpf(float(z))
            factor=1+b*(mp.mpf(float(eta))-mp.mpf('-0.00008405268308670538'))/2
            if factor<=0:
                raise ValueError('Nonpositive quadratic mass-law factor')
            x=r*factor*factor
            xphi=b*r*factor
            xphiphi=b*b*r/2
            value=self.action.action_sources(x,zz,r,dps=self.dps,N=self.terms)
            q=value['q'];u=x/zz;a=mp.mpf(3)/2;c0,c1,c2,ak=self.action.coefficients(self.dps,self.terms)
            B0=mp.zeta(-3,a)-mp.zeta(-1,a)/4+q/8+q*q/4
            B0q=mp.mpf(1)/8+q/2
            B1=c0+c1*q+c2*q*q
            B1q=c1+2*c2*q
            B1qq=2*c2
            for k,A in enumerate(ak,3):
                B1+=(-q)**k*A/k
                B1q+=(-1)**k*q**(k-1)*A
                B1qq+=(-1)**k*(k-1)*q**(k-2)*A
            P=B1-B0*mp.log(zz)
            Pq=B1q-B0q*mp.log(zz)
            Pqq=B1qq-mp.log(zz)/2
            C0=1/(16*mp.pi**2)
            Wxx=-C0*Pqq-C0*mp.log(r)/2
            Wxz=-C0*(Pq-u*Pqq-B0q)+C0*mp.log(r)
            Wzz=-C0*(2*P-2*u*Pq-3*B0+u*u*Pqq+2*u*B0q)-C0*mp.mpf(29)/15*mp.log(r)
            rhox=value['Wx']-zz*Wxz/2
            rhoz=value['Wy']/2-zz*Wzz/2
            Qx=2*Wxx;Qz=2*Wxz
            j=xphi*value['Q']/2
            theta=abs(q)/(a*a)
            common=a**3*(1+a/(2*self.terms-2))*theta**(self.terms+1)/(1-theta)
            eB=common/(self.terms+1)
            eD=common/abs(q) if q else mp.mpf(0)
            eDD=common/(q*q)*(self.terms+theta/(1-theta)) if q else mp.mpf(0)
            eWx=C0*zz*eD;eWy=C0*zz*(2*eB+abs(u)*eD)
            eWxz=C0*(eD+abs(u)*eDD)
            eWzz=C0*(2*eB+2*abs(u)*eD+u*u*eDD)
            eQx=2*C0*eDD;eQz=2*eWxz
            erho_x=eWx+zz*eWxz/2;erho_z=eWy/2+zz*eWzz/2
            derivative_bounds={'rho_x':erho_x,'rho_z':erho_z,'Q_x':eQx,'Q_z':eQz,
                'rho_eta':abs(xphi)*erho_x,
                'j_eta':(abs(xphiphi)*value['tail_bound']['Q']+xphi*xphi*eQx)/2,
                'j_z':abs(xphi)*eQz/2}
            result={k:float(value[k]) for k in ['W','rho','Q','Wx','Wy']}
            result.update(x=float(x),x_phi=float(xphi),x_phiphi=float(xphiphi),mass_factor=float(factor),
                          u=float(u),j=float(j),rho_eta=float(xphi*rhox),
                          rho_z=float(rhoz),j_eta=float((xphiphi*value['Q']+xphi*xphi*Qx)/2),
                          j_z=float(xphi*Qz/2),rho_x=float(rhox),Q_x=float(Qx),Q_z=float(Qz),
                          derivative_tail_bound={k:float(v) for k,v in derivative_bounds.items()},
                          tail_bound={k:float(v) for k,v in value['tail_bound'].items()})
        self.source_cache[key]=result
        return result

    def integrate(self, u, resolution):
        ell,yb=map(float,u)
        if not (-28<ell<-18 and 25<yb<35):
            raise ValueError('Shooting coordinates outside registered search box')
        eh=-10.**ell
        start=resolution['y_start']
        sol=solve_ivp(self.bulk.rhs,(start,yb),self.bulk.cone(eh,start,True),
                      method=resolution['method'],rtol=resolution['rtol'],
                      atol=[1e-14]+[1e-100]*6,max_step=resolution['max_step'],dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        sampled=np.unique(np.concatenate((sol.t,np.linspace(start,yb,1001))))
        profile=sol.sol(sampled)
        Up=self.bulk.U(profile[1]);wp=profile[2]
        scales={'sample_count':len(sampled),
                'max_sampled_sqrt_abs_kT':float(np.sqrt(np.max(np.abs(Up/6-wp*wp/12)))),
                'max_sampled_sqrt_abs_kN':float(np.sqrt(np.max(np.abs(Up/6+wp*wp/4))))}
        return sol.y[:,-1],sol.nfev,scales

    def endpoint(self, u, b, gamma, resolution, mutation=None):
        state,nfev,bulk_scales=self.integrate(u,resolution)
        R,e,w,A,f,g,I=state
        z=1/R**2
        if not (.75<=z/Z0<=1.25 and abs(e-ETA0)<=.1*abs(ETA0)):
            raise ValueError('Endpoint outside registered response domain')
        source=self.source(e,z,b)
        if not (0<source['u']<4.5):
            raise ValueError('Positive-gap convergent source domain violated')
        evaluated_source=dict(source)
        if mutation=='freeze_sources':
            frozen=self.source(ETA0,Z0,b)
            source=dict(source, rho=frozen['rho'], j=frozen['j'],rho_eta=0.,rho_z=0.,j_eta=0.,j_z=0.)
        if mutation=='omit_metric_source':
            source=dict(source,rho=0.,rho_eta=0.,rho_z=0.)
        if mutation=='flip_current':
            source=dict(source,j=-source['j'],j_eta=-source['j_eta'],j_z=-source['j_z'])
        W=1/3+e*e+e**3/3
        p=2*e+e*e
        sigma=self.bulk.sigma(e)
        sp=self.bulk.sig1(e)
        s=gamma*source['rho'];t=gamma*source['j']
        tau=DELTA*(1+DETUNING_C*(1+e))+s
        # Exact constraint after analytic cancellation of the superpotential square.
        C=z+(w-p)*(w+p)/12-W*tau/9-tau*tau/36
        B=w+sp/2+t/2
        q=math.sqrt(z+w*w/12-self.bulk.U(e)/6)
        phi_a=np.array([f,w])
        z_a=np.array([-2*A/R**3,-2*z*q])
        q_a=np.array([-2*I/R**4-w*f/3,-z-w*w/3])
        wp=self.bulk.U(e,1)-4*q*w
        w_a=np.array([g,wp])
        s_a=gamma*(source['rho_eta']*phi_a+source['rho_z']*z_a)
        t_a=gamma*(source['j_eta']*phi_a+source['j_z']*z_a)
        M1=q_a-sp*phi_a/6-s_a/6
        M2=w_a+self.bulk.sig2(e)*phi_a/2+t_a/2
        E1=C/(q+(sigma+s)/6)
        # Stable integrated q_ell identity, with the exact off-shell product rule.
        C_a=(q+(sigma+s)/6)*M1+E1*(q_a+(sp*phi_a+s_a)/6)
        C_direct=z_a+w*w_a/6-self.bulk.U(e,1)*phi_a/6-(sigma+s)*(sp*phi_a+s_a)/18
        J=np.array([C_a,M2])
        original=np.array([E1,B])
        energy=max(math.sqrt(z),math.sqrt(evaluated_source['x']),abs(q),
                   bulk_scales['max_sampled_sqrt_abs_kT'],bulk_scales['max_sampled_sqrt_abs_kN'])
        hierarchy={'Ehat':energy,'energy_over_M5':energy*gamma**(1/3),
                   'tenfold_screen_pass':bool(energy*gamma**(1/3)<=.1),
                   'status':'classical_control' if gamma==0 else 'declared_scale_screen',**bulk_scales}
        return dict(coordinates=list(map(float,u)),eh=float(-10.**u[0]),yb=float(u[1]),
            R=float(R),eta=float(e),phi=float(1+e),w=float(w),H2=float(z),q=float(q),
            sigma=float(sigma),sigma_phi=float(sp),source=evaluated_source,
            applied_source=source,s=float(s),t=float(t),C=float(C),B=float(B),
            residual=original.tolist(),scaled_residual=(np.array([C,B])/SCALES).tolist(),
            jacobian=J.tolist(),jacobian_constraint_discrepancy=(C_a-C_direct).tolist(),
            jacobian_condition=float(np.linalg.cond(J/SCALES[:,None])),
            jacobian_determinant=float(np.linalg.det(J)),
            source_jacobian=np.array([s_a,t_a]).tolist(),output_jacobian=np.array([phi_a,z_a]).tolist(),
            state=list(map(float,state)),nfev=nfev,hierarchy=hierarchy,positive_tension=bool(sigma+s>0),
            relative_H2_change=float(z/Z0-1),relative_eta_change=float((e-ETA0)/abs(ETA0)))

    def solve(self, seed, b, gamma, resolution, mutation=None):
        u=np.array(seed,dtype=float)
        history=[]
        evaluations=0
        began=time.monotonic()
        rejected_trials=[]
        for iteration in range(12):
            value=self.endpoint(u,b,gamma,resolution,mutation)
            evaluations+=1
            residual=np.array([value['C'],value['B']])
            norm=float(np.linalg.norm(residual/SCALES,np.inf))
            history.append({'iteration':iteration,'coordinates':u.tolist(),'C':value['C'],'B':value['B'],'scaled_norm':norm})
            if abs(value['C'])<=2e-18 and abs(value['B'])<=2e-14:
                value.update(b=b,gamma=gamma,resolution=resolution['name'],mutation=mutation,
                             iterations=history,rejected_trials=rejected_trials,
                             evaluations=evaluations,seconds=time.monotonic()-began)
                return value
            step=np.linalg.solve(np.array(value['jacobian']),-residual)
            limit=max(1.,abs(step[0])/.3,abs(step[1])/.2)
            step/=limit
            accepted=False
            for halvings in range(12):
                trial=u+step*2.**(-halvings)
                try:
                    candidate=self.endpoint(trial,b,gamma,resolution,mutation)
                    evaluations+=1
                    trialnorm=np.linalg.norm(np.array([candidate['C'],candidate['B']])/SCALES,np.inf)
                except (ValueError,RuntimeError) as exc:
                    rejected_trials.append({'iteration':iteration,'halvings':halvings,
                                            'coordinates':trial.tolist(),'error':str(exc)})
                    continue
                if trialnorm<norm:
                    u=trial;accepted=True;break
                rejected_trials.append({'iteration':iteration,'halvings':halvings,
                                        'coordinates':trial.tolist(),'scaled_norm':float(trialnorm),
                                        'reason':'no decrease'})
            if not accepted:
                raise RuntimeError('Newton backtracking failed: '+json.dumps({'iterations':history,'rejected_trials':rejected_trials}))
        raise RuntimeError('Newton iteration cap reached: '+json.dumps({'iterations':history,'rejected_trials':rejected_trials}))

def run(repo, output, public_freeze):
    protocol=verify_inputs(repo)
    if not public_freeze:
        raise RuntimeError('An explicit public freeze commit/link is required')
    experiment=protocol['experiment']
    output.mkdir(parents=True,exist_ok=False)
    record={'started_at_utc':utc(),'public_freeze':public_freeze,
        'registration_sha256':sha(ROOT/'PREREGISTRATION.json'),
        'source_commit':protocol['source_commit'],
        'observed_checkout_head':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
        'producer_sha256':sha(__file__),'python':platform.python_version(),
        'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,
        'runs':[],'negative_controls':[],'failures':[]}
    write_json(output/'started.json',{k:v for k,v in record.items() if k not in ['runs','negative_controls','failures']})
    model=Model(repo)
    archived_u=np.array([math.log10(-EH0),Y0])
    response=np.array(model.archive['coordinate_response'])
    obs_response=np.array(model.archive['observable_response'])
    baseline_sources={b:model.source(ETA0,Z0,b) for b in experiment['b_values']}
    record['baseline_sources']={str(b):v for b,v in baseline_sources.items()}
    fine_by_case={}
    for resolution in experiment['resolutions']:
        for b in experiment['b_values']:
            last=archived_u.copy();last_gamma=0.
            grid=experiment['gamma_values'] if resolution['coverage']=='all' else [0.,experiment['gamma_values'][-1]]
            for gamma in grid:
                case={'resolution':resolution['name'],'b':b,'gamma':gamma}
                src=baseline_sources[b]
                if (b,gamma) in fine_by_case:
                    seed=np.array(fine_by_case[(b,gamma)]['coordinates'])
                else:
                    seed=last+response@np.array([src['rho'],src['j']])*(gamma-last_gamma)
                try:
                    result=model.solve(seed,b,gamma,resolution)
                    result['linear_prediction_coordinates']=(archived_u+gamma*response@np.array([src['rho'],src['j']])).tolist()
                    result['linear_prediction_delta_observables']=(gamma*obs_response@np.array([src['rho'],src['j']])).tolist()
                    record['runs'].append(result)
                    last=np.array(result['coordinates']);last_gamma=gamma
                    if resolution['name']=='refined':fine_by_case[(b,gamma)]=result
                    print('SOLVED',case,'C=',result['C'],'B=',result['B'],flush=True)
                except Exception:
                    record['failures'].append({'case':case,'seed':seed.tolist(),'traceback':traceback.format_exc()})
                    print('FAILED',case,flush=True)
                write_json(output/'results.json',record)
    resolution=next(r for r in experiment['resolutions'] if r['name']=='refined')
    gamma=experiment['gamma_values'][-1]
    for mutation in experiment['negative_controls']:
        case={'b':1,'gamma':gamma,'mutation':mutation}
        try:
            seed=fine_by_case[(1,gamma)]['coordinates']
            wrong=model.solve(seed,1,gamma,resolution,mutation)
            correct=model.endpoint(wrong['coordinates'],1,gamma,resolution)
            wrong['correct_equation_residual']=[correct['C'],correct['B']]
            record['negative_controls'].append(wrong)
            print('CONTROL',mutation,wrong['correct_equation_residual'],flush=True)
        except Exception:
            record['failures'].append({'case':case,'traceback':traceback.format_exc()})
        write_json(output/'results.json',record)
    record['completed_at_utc']=utc()
    record['execution_status']='COMPLETE' if not record['failures'] else 'COMPLETE_WITH_FAILURES'
    write_json(output/'results.json',record)
    return 0 if not record['failures'] else 1

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo',type=Path,default=Path('/workspace/HDblast'))
    ap.add_argument('--freeze',action='store_true')
    ap.add_argument('--output',type=Path)
    ap.add_argument('--public-freeze',help='Published preregistration commit or permanent URL')
    args=ap.parse_args()
    if args.freeze:
        print(json.dumps(freeze(args.repo),indent=2));return 0
    if args.output is None:
        ap.error('--output is required')
    return run(args.repo,args.output,args.public_freeze)

if __name__=='__main__':
    sys.exit(main())
