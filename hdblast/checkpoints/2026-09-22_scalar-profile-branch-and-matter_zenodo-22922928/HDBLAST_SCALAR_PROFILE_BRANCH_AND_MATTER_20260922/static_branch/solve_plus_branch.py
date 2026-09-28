"""Floating regular-cone AdS-plus static branch satisfying BOTH junctions.

Uses eta=phi-1 throughout the very small cone displacement. This is an
existence candidate, not an interval certificate or an attractor proof.
"""
from pathlib import Path
import json, math, hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.special import eval_gegenbauer

C=2/1.0357712571566784-4/3

def pot(e):
    w=1/3+e*e+e**3/3;wp=e*(2+e);p=1+e
    u=.5*wp*wp-(2/3)*w*w
    up=wp*(2*p-(4/3)*w)
    upp=4*p*p+2*wp-(4/3)*(wp*wp+2*p*w)
    return u,up,upp

def solve(delta,rtol=2e-12,y0=1e-4):
    sig0=2/3+delta*(1+C)
    hmetric2=sig0*sig0/36-1/81
    yguess=9*np.arcsinh(1/(9*np.sqrt(hmetric2)))
    x=np.cosh(yguess/9)
    f=eval_gegenbauer(14,2,x)/eval_gegenbauer(14,2,1.)
    g=(1/9)*np.sinh(yguess/9)*4*eval_gegenbauer(13,3,x)/eval_gegenbauer(14,2,x)
    eguess=-delta*C/(2*(g+2))
    def integrate(params,dense=False):
        eh=-10**params[0];yb=params[1];u,up,upp=pot(eh)
        aa=-u/36;bb=u*u/4320-up*up/750;cc=up/10;dd=up*(upp/280+u/630)
        ini=[y0+aa*y0**3+bb*y0**5,eh+cc*y0*y0+dd*y0**4,2*cc*y0+4*dd*y0**3]
        def rhs(y,v):
            rho,eta,py=v;u,up,_=pot(eta)
            rad=1+rho*rho*(py*py/12-u/6)
            if rad<=0:raise ValueError('Nonpositive radial first-integral branch')
            ry=np.sqrt(rad)
            return [ry,py,up-4*ry*py/rho]
        sol=solve_ivp(rhs,(y0,yb),ini,method='DOP853',rtol=rtol,atol=[1e-14,1e-42,1e-42],
                      max_step=.1,dense_output=dense)
        if not sol.success:raise RuntimeError(sol.message)
        return sol
    def residual(params):
        rho,e,py=integrate(params).y[:,-1];u,_,_=pot(e)
        ry=np.sqrt(1+rho*rho*(py*py/12-u/6));w=1/3+e*e+e**3/3
        sig=2*w+delta*(1+C*(1+e));sig1=2*e*(2+e)+delta*C
        return np.array([ry/rho-sig/6,py+sig1/2])/delta
    guess=[math.log10(abs(eguess/f)),yguess]
    fit=root(residual,guess,tol=1e-9,options={'eps':1e-8})
    res=residual(fit.x)*delta
    if np.max(np.abs(res))>1e-10:raise RuntimeError('Junction root did not converge')
    sol=integrate(fit.x,True);rho,e,py=sol.y[:,-1];u,_,_=pot(e)
    w=1/3+e*e+e**3/3;sig=2*w+delta*(1+C*(1+e));sig1=2*e*(2+e)+delta*C
    h2=1/rho**2;static_h2=sig*sig/36-sig1*sig1/48+u/6
    series2=delta*(1+C)/27+delta**2*((1+C)**2/36-C*C/384)
    ys=np.linspace(y0,fit.x[1],4097);arr=sol.sol(ys)
    return dict(delta=delta,rtol=rtol,y0=y0,root_success=bool(fit.success),eta_h=-10**fit.x[0],
                y_b=float(fit.x[1]),rho_b=float(rho),eta_b=float(e),phi_b=float(1+e),phi_y_b=float(py),
                junction_residual=res.tolist(),H2=float(h2),junction_H2=float(static_h2),
                H2_identity_residual=float(h2-static_h2),metric_only_H2=float(hmetric2),
                H2_second_order_expansion=float(series2),H2_expansion_remainder=float(h2-series2),
                eta_first_order=-9*delta*C/64,eta_finite_curvature_linear=float(eguess),
                constant_phi_scalar_junction_residual=delta*C/2,
                min_phi=float(1+arr[1].min()),max_abs_phi_y=float(np.max(np.abs(arr[2]))),
                scope='Floating regular-cone static branch candidate; no dynamical attraction, stability or interval existence claim'),np.vstack([ys,arr])

if __name__=='__main__':
    out=Path(__file__).resolve().parent;rows=[];profiles={}
    for d in [.0003,.001,.003,.01,.03,.1]:
        row,p=solve(d);ref,p2=solve(d,rtol=8e-14,y0=5e-5);row['refined']=ref
        rows.append(row);profiles[str(d)]=p2;print(json.dumps(row),flush=True)
    result=dict(status='FLOATING_STATIC_PLUS_BRANCH',c=C,rows=rows,
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'PLUS_BRANCH_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(out/'PLUS_BRANCH_PROFILES.npz',**profiles)
