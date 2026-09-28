#!/usr/bin/env python3
"""Recompute nonlinear constraint term balances from saved states; no PDE runs.

NumPy only. Uses explicit finite differences, separately constructed polynomial
ghosts, and analytic reference derivatives. Reports original and actual-term
normalizations side by side, never silently replacing the producer's norm.
"""
from pathlib import Path
import hashlib,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def W(p):return 1-p+p**3/3
def U(p):return .5*(p*p-1)**2-(2/3)*W(p)**2
C=2/1.0357712571566784-4/3

def analyze(path):
    meta=json.loads(path.read_text());archive=path.with_suffix('.npz')
    with np.load(archive,allow_pickle=False) as f:
        z=f['z'].copy();rho,phi,hc,phz=f['background']; snapshots={k:f[k].copy() for k in f.files if k.startswith('state_')}
        stored={k:f[k].copy() for k in f.files if k.startswith('constraints_')}
    n=len(z);h=meta['hmin'];stretch=meta.get('stretch',.05);zp=h*(1-z/stretch);zpp=-h*zp/stretch
    rb=rho[-1];pb=phi[-1];sig0=2*W(pb)+.001*(1+C*pb);sig10=2*(pb*pb-1)+.001*C;sig20=4*pb
    mat=np.zeros((6,6));mat[:5]=(-np.arange(5.))[:,None]**np.arange(6)[None,:];mat[5,1]=1.
    def bc(b,f,bt,ft):
        eb=np.exp(b);ds=sig10*f+2*pb*f*f+(2/3)*f**3;ds1=sig20*f+2*f*f
        ga=rb*(sig0*np.expm1(b)+eb*ds)/6;gf=-rb*(sig10*np.expm1(b)+eb*ds1)/2
        gat=rb*eb*((sig0+ds)*bt+(sig10+ds1)*ft)/6
        return ga,gf,gat
    def derivatives(v,g):
        coeff=np.linalg.solve(mat,np.r_[v[-1:-6:-1],g*zp[-1]])
        ghost=np.polynomial.polynomial.polyval(np.arange(1.,3.),coeff)
        ext=np.r_[0.,0.,v,ghost]
        fx=(ext[:-4]-8*ext[1:-3]+8*ext[3:-1]-ext[4:])/12
        fxx=(-ext[:-4]+16*ext[1:-3]-30*ext[2:-2]+16*ext[3:-1]-ext[4:])/12
        return fx/zp,fxx/zp**2-fx*zpp/zp**3
    mask=(z>-.8*meta['L'])&(np.arange(n)>5);ids=np.flatnonzero(mask);core=z>-.2
    scale=1+6*hc*hc+phz*phz
    rows=[];initial_max=None;initial_fields=None;previous=None
    for key in sorted(snapshots,key=lambda x:int(x.split('_')[1])):
        step=int(key.split('_')[1]);t=step*meta['dt'];a,b,f,pa,pbt,pf=snapshots[key]
        ga,gf,gat=bc(b[-1],f[-1],pbt[-1],pf[-1]);az,azz=derivatives(a,ga);bz,_=derivatives(b,ga);fz,_=derivatives(f,gf);paz,_=derivatives(pa,gat)
        paz[-1]=np.dot(np.array([25/12,-4,3,-4/3,1/4]),pa[-1:-6:-1])/zp[-1]
        At=1+pa;Az=hc+az;Bt=pbt;Bz=hc+bz;Pz=phz+fz;Azz=hc*hc-1-phz*phz/3+azz
        termH=np.array([-2*rho*rho*np.exp(2*b)*U(phi+f),6*At*At,6*At*Bt,-12*Az*Az,6*Az*Bz,-6*Azz,-pf*pf,-Pz*Pz])
        termM=np.array([-3*paz,-3*At*Az,3*At*Bz,3*Az*Bt,-pf*Pz])
        H=termH.sum(axis=0);M=termM.sum(axis=0);sumH=np.abs(termH).sum(axis=0);sumM=np.abs(termM).sum(axis=0)
        nrH=np.abs(H)/np.maximum(sumH,1e-300);nrM=np.abs(M)/np.maximum(sumM,1e-300)
        weight=(rho/rb)**3*np.exp(3*(t+a));Cp=weight*(H+2*M);Cm=weight*(H-2*M)
        maxC=max(float(np.max(abs(Cp[mask]))),float(np.max(abs(Cm[mask]))))
        if initial_max is None:initial_max=maxC;initial_fields=(Cp.copy(),Cm.copy())
        def point(j):
            return dict(z=float(z[j]),local_grid_scale=float(zp[j]),H=float(H[j]),M=float(M[j]),
                original_H_norm=float(abs(H[j])/scale[j]),original_M_norm=float(abs(M[j])/scale[j]),
                actual_H_term_sum=float(sumH[j]),actual_M_term_sum=float(sumM[j]),
                actual_H_relative_residual=float(nrH[j]),actual_M_relative_residual=float(nrM[j]),
                H_terms=termH[:,j].tolist(),M_terms=termM[:,j].tolist(),
                weighted_outgoing=float(Cp[j]),weighted_incoming=float(Cm[j]),weight=float(weight[j]))
        jH=ids[np.argmax(abs(H[mask])/scale[mask])];jAH=ids[np.argmax(nrH[mask])];jC=ids[np.argmax(np.maximum(abs(Cp[mask]),abs(Cm[mask])))]
        ray=(z>=-t-.05)&(z<=min(0.,-t+.025))&mask
        ray_ids=np.flatnonzero(ray); raydiag=None
        if len(ray_ids):
            jr=ray_ids[np.argmax(abs(Cp[ray]))]
            backward=z[jr]+t
            expected=float(np.interp(backward,z,initial_fields[0])) if z[0]<=backward<=0 else None
            raydiag=dict(outgoing_peak=point(jr),initial_characteristic_foot=float(backward),
                         transported_initial_Cplus=expected,
                         departure_from_initial_transport=None if expected is None else float(Cp[jr]-expected))
            if previous is not None:
                tp,oldCp,oldCm=previous;foot=z[jr]+t-tp
                oldval=float(np.interp(foot,z,oldCp)) if z[0]<=foot<=0 else None
                raydiag['previous_snapshot_time']=tp
                raydiag['previous_characteristic_foot']=float(foot)
                raydiag['advected_previous_Cplus']=oldval
                raydiag['departure_from_previous_transport']=None if oldval is None else float(Cp[jr]-oldval)
        # Independently reconstructed H/M should match the saved deviation form,
        # allowing cancellation noise from our full-term evaluation.
        saved=stored.get('constraints_'+str(step));repro=None
        if saved is not None:repro=dict(H_max_abs_difference=float(np.max(abs(H[mask]-saved[0,mask]))),M_max_abs_difference=float(np.max(abs(M[mask]-saved[1,mask]))))
        rows.append(dict(time=t,original_H_norm_max=float(np.max(abs(H[mask])/scale[mask])),
                         original_M_norm_max=float(np.max(abs(M[mask])/scale[mask])),
                         actual_H_term_relative_max=float(np.max(nrH[mask])),actual_M_term_relative_max=float(np.max(nrM[mask])),
                         core_actual_H_relative_max=float(np.max(nrH[core])),core_actual_M_relative_max=float(np.max(nrM[core])),
                         peak_original_H_norm=point(jH),peak_actual_H_relative=point(jAH),
                         peak_weighted_characteristic=point(jC),weighted_C_max=maxC,
                         weighted_C_max_over_initial=maxC/max(initial_max,1e-300),outgoing_ray=raydiag,
                         reproduction_difference=repro))
        previous=(t,Cp.copy(),Cm.copy())
    return dict(run=path.stem,epsilon=meta['epsilon'],hmin=h,stretch=stretch,L=meta['L'],nodes=n,
                source_json_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),source_npz_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),snapshots=rows)

runs=[]
for p in sorted((ROOT/'runs').glob('balanced_*.json')):
    if p.with_suffix('.npz').exists():runs.append(analyze(p))
out=dict(status='READ_ONLY_SAVED_STATE_RESIDUAL_ANALYSIS',scope='No evolution. Actual-term denominator is an additional cancellation ratio, not a replacement for original diagnostics or an error certificate.',
         H_term_order=['-2exp(2B)U','6At^2','6AtBt','-12Az^2','6AzBz','-6Azz','-phi_t^2','-phi_z^2'],
         M_term_order=['-3Atz','-3AtAz','3AtBz','3AzBt','-phi_t phi_z'],runs=runs)
Path(__file__).with_name('EVOLUTION_RESIDUAL_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
for run in runs:
    last=run['snapshots'][-1]
    print(json.dumps({'run':run['run'],'time':last['time'],'H_background_norm':last['original_H_norm_max'],
        'H_actual_terms':last['actual_H_term_relative_max'],'M_actual_terms':last['actual_M_term_relative_max'],
        'core_H_actual_terms':last['core_actual_H_relative_max'],'core_M_actual_terms':last['core_actual_M_relative_max'],
        'Hpeak_z':last['peak_original_H_norm']['z'],'weighted_growth':last['weighted_C_max_over_initial'],
        'weighted_peak_z':last['peak_weighted_characteristic']['z'],'reproduction':last['reproduction_difference']}),flush=True)
