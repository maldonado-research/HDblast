#!/usr/bin/env python3
"""Bounded independent checks of the continuation solver; no evolution runs.

Requires NumPy and SciPy. Reads existing saved spectral modes and checks ghost
polynomials, the analytic map, Jacobian, and mode-constraint amplitude scaling.
"""
from pathlib import Path
import importlib.util,json,hashlib
import numpy as np
root=Path(__file__).resolve().parents[1]
solver_path=root/'solver/registered_solver.py'
spec=importlib.util.spec_from_file_location('new_solver_audit_target',solver_path)
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
checks=[]
def check(name,ok,**extra):
    checks.append(dict(name=name,pass_check=bool(ok),**extra))
    if not ok:raise RuntimeError(name)

gh=R.hermite_ghosts()
for degree in range(6):
    vals=(-np.arange(5.))**degree
    result=gh@np.r_[vals,1. if degree==1 else 0.]
    error=float(np.max(np.abs(result-np.arange(1.,4.)**degree)))
    check('Hermite reproduces degree '+str(degree),error<3e-10,error=error)

z,zp,zpp=R.grid(.0002,6,.05)
alpha=.0002/.05
check('analytic grid orientation',bool(np.all(np.diff(z)>0) and np.all(zp>0) and np.all(zpp<0)))
check('analytic grid second derivative sign',np.max(np.abs(zpp+alpha*zp))==0.)
D1,D2,g1,g2,_=R.matrices(zp,zpp)
for k in [0,1,2]:
    q=z**k; slope=1. if k==1 else 0.
    d1=D1@q+g1*slope; d2=D2@q+g2*slope
    e1=np.zeros_like(z) if k==0 else k*z**(k-1)
    e2=np.zeros_like(z) if k<2 else k*(k-1)*z**(k-2)
    # Exclude the clamped left boundary but include the physical shell.
    rel1=float(np.max(np.abs(d1[4:]-e1[4:])/(1+np.abs(e1[4:]))))
    rel2=float(np.max(np.abs(d2[4:]-e2[4:])/(1+np.abs(e2[4:]))))
    check('mapped derivative of z^'+str(k),rel1<1e-7 and rel2<1e-6,relative_error_D1=rel1,relative_error_D2=rel2)

def restored_solver(npz,meta):
    obj=R.Solver.__new__(R.Solver)
    obj.z,obj.zp,obj.zpp=R.grid(meta['hmin'],meta['L'],meta['stretch']);obj.n=len(obj.z)
    obj.D1,obj.D2,obj.g1,obj.g2,obj.ghost=R.matrices(obj.zp,obj.zpp)
    obj.rho,obj.phi,obj.hc,obj.phz=npz['background']
    obj.rho2=obj.rho**2;obj.rb=obj.rho[-1];obj.pb=obj.phi[-1]
    obj.tdet=meta['background']['tdet'];obj.ko=0.;obj.hmin=meta['hmin'];obj.L=meta['L'];obj.stretch=meta['stretch']
    obj.pot=R.derivatives(obj.phi);obj.s0=2*R.W(obj.pb)+obj.tdet*(1+R.C*obj.pb)
    obj.s10=2*(obj.pb**2-1)+obj.tdet*R.C;obj.s20=4*obj.pb
    if np.max(np.abs(obj.z-npz['z']))>1e-12:raise RuntimeError('grid mismatch')
    return obj

rows=[]
for name in ['spectrum_reg_h4','spectrum_reg_h2','spectrum_reg_h1','spectrum_reg_h05']:
    mp=root/'runs'/(name+'.json'); npz=root/'runs'/(name+'.npz')
    if not mp.exists() or not npz.exists():continue
    meta=json.loads(mp.read_text())
    with np.load(npz,allow_pickle=False) as saved:
        obj=restored_solver(saved,meta);mode=saved['mode'].copy()
    check(name+' mode has negligible imaginary eigenvalue',abs(meta['eigenvalues'][0]['imag'])<1e-12)
    amp=1e-8
    _,hp,mp_,cp=obj.diagnostics(amp*mode,0.)
    _,hn,mn,cn=obj.diagnostics(-amp*mode,0.)
    hlin=(hp-hn)/(2*amp);mlin=(mp_-mn)/(2*amp);clin=(cp-cn)/(2*amp)
    heven=(hp+hn)/(2*amp**2);meven=(mp_+mn)/(2*amp**2)
    mask=(obj.z>-.8*obj.L)&(np.arange(obj.n)>5)
    ids=np.flatnonzero(mask);j=ids[np.argmax(np.max(np.abs(clin[:,mask]),axis=0))]
    jb=ids[np.argmax(np.abs(hlin[mask])+np.abs(mlin[mask]))]
    # Independent linear H and M via centered/one-sided derivatives in xi,
    # not inserting junction derivatives. 4th-order in the interior; one-sided
    # polynomial stencils at the two rightmost nodes.
    def independent_d1(v):
        d=np.zeros_like(v);d[2:-2]=(v[:-4]-8*v[1:-3]+8*v[3:-1]-v[4:])/12
        for j in [-2,-1]:
            at=obj.n+j;idx=np.arange(obj.n-5,obj.n);xx=idx-at
            V=xx[:,None]**np.arange(5)[None,:]
            weights=np.linalg.solve(V.T,np.eye(5)[1]);d[at]=weights@v[idx]
        return d/obj.zp
    def independent_d2(v):
        d=np.zeros_like(v);dx=np.zeros_like(v)
        d[2:-2]=(-v[:-4]+16*v[1:-3]-30*v[2:-2]+16*v[3:-1]-v[4:])/12
        dx[2:-2]=(v[:-4]-8*v[1:-3]+8*v[3:-1]-v[4:])/12
        for j in [-2,-1]:
            at=obj.n+j;idx=np.arange(obj.n-6,obj.n);xx=idx-at
            V=xx[:,None]**np.arange(6)[None,:]
            w2=np.linalg.solve(V.T,2*np.eye(6)[2]);w1=np.linalg.solve(V.T,np.eye(6)[1]);d[at]=w2@v[idx];dx[at]=w1@v[idx]
        return d/obj.zp**2-dx*obj.zpp/obj.zp**3
    a,b,f,pa,pb,pf=mode;az=independent_d1(a);bz=independent_d1(b);fz=independent_d1(f);paz=independent_d1(pa)
    independent_H=-2*obj.rho2*(2*b*obj.pot[0]+f*obj.pot[1])+12*pa+6*pb-18*obj.hc*az+6*obj.hc*bz-6*independent_d2(a)-2*obj.phz*fz
    independent_M=-3*paz-3*(az-bz)+3*obj.hc*pb-pf*obj.phz
    one=np.array([25/12,-4,3,-4/3,1/4])/obj.zp[-1]
    slope_defects=[float(one@mode[i,-1:-6:-1]-target) for i,target in zip([0,1,2],[obj.rb*(obj.s0*b[-1]+obj.s10*f[-1])/6]*2+[-obj.rb*(obj.s10*b[-1]+obj.s20*f[-1])/2])]
    ev=meta['eigenvalues'][0]['real']; projected=mode.copy();projected[3:]=ev*mode[:3]
    _,hhp,mmp,ccp=obj.diagnostics(amp*projected,0.)
    _,hhn,mmn,ccn=obj.diagnostics(-amp*projected,0.)
    linear_M_projected=(mmp-mmn)/(2*amp)
    ga_quadratic=obj.rb*(.5*obj.s0*b[-1]**2+obj.s10*b[-1]*f[-1]+2*obj.pb*f[-1]**2)/6
    leading_H_corner=-6*obj.g2[-1]*ga_quadratic
    rows.append(dict(name=name,hmin=meta['hmin'],eigenvalue=ev,
        mode_max_absolute_by_field=np.max(np.abs(mode),axis=1).tolist(),
        per_unit_shell_scalar_amplitude=dict(weighted_characteristic_max=float(np.max(np.abs(clin[:,mask]))),
            weighted_peak_z=float(obj.z[j]),weighted_peak_index=int(j),
            H_at_weighted_peak=float(hlin[j]),M_at_weighted_peak=float(mlin[j]),
            maximum_unweighted_H_plus_M_position=float(obj.z[jb]),
            independent_shell_H=float(independent_H[-1]),independent_shell_M=float(independent_M[-1]),
            independent_shell_field_slope_defects=slope_defects),
        velocity_eigenvector_defect_max=float(np.max(np.abs(mode[3:]-ev*mode[:3]))),
        projected_velocities_diagnostic_only=dict(shell_linear_M=float(linear_M_projected[-1]),shell_linear_M_before=float(mlin[-1])),
        second_order_even_residual_max={'H':float(np.max(np.abs(heven[mask]))),'M':float(np.max(np.abs(meven[mask])))},
        leading_second_order_H_corner_from_nonlinear_boundary=float(leading_H_corner)))
    # Independent central difference of nonlinear RHS vs assembled sparse Jacobian,
    # complementing the producer's single complex-step directional check.
    rng=np.random.default_rng(199);v=rng.normal(size=mode.shape);v[:,:2]=0
    eps=1e-6;fd=(obj.rhs(eps*v)-obj.rhs(-eps*v))/(2*eps);exact=(obj.linear_matrix()@v.ravel()).reshape(v.shape)
    err=float(np.max(np.abs(fd-exact))/(1+np.max(np.abs(exact))))
    check(name+' real directional Jacobian',err<1e-9,relative_error=err)

result=dict(status='PASS',scope='No evolution or spectrum recomputed. New polynomial, map, real Jacobian and saved-mode constraint diagnostics.',
            solver_sha256=hashlib.sha256(solver_path.read_bytes()).hexdigest(),checks=checks,mode_diagnostics=rows)
Path(__file__).with_name('NEW_SOLVER_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'modes':rows},indent=2))
