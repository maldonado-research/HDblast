"""Complex DOP853 smooth-step controls and exact UV cutoff diagnostics."""
import json, hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
def logsinh(x):
    x=np.asarray(x)
    return x+np.log(-np.expm1(-2*x))-np.log(2.)
def spectrum(k,eps):
    wm=np.sqrt(k*k+1);wp=np.sqrt(k*k+4)
    return np.exp(2*logsinh(np.pi*eps*(3/(wp+wm))/2)-logsinh(np.pi*eps*wm)-logsinh(np.pi*eps*wp))
def run():
    cases=[]
    for eps in [.1,.3,1.]:
        for k in [.5,2.,8.]:
            wm=np.sqrt(k*k+1);wp=np.sqrt(k*k+4);lo=-18*eps;hi=18*eps
            u=np.exp(-1j*wm*lo)/np.sqrt(2*wm);v=-1j*wm*u
            def rhs(t,z):return np.array([z[1],-(k*k+2.5+1.5*np.tanh(t/eps))*z[0]])
            sol=solve_ivp(rhs,(lo,hi),[u,v],method="DOP853",rtol=2e-11,atol=2e-13,
                          max_step=min(eps/20,.2/wp))
            if not sol.success:raise RuntimeError(sol.message)
            u,v=sol.y
            beta=(wp*u[-1]-1j*v[-1])/np.sqrt(2*wp)*np.exp(-1j*wp*hi)
            num=float(abs(beta)**2);ex=float(spectrum(k,eps))
            err=abs(num-ex);rel=err/ex
            wr=float(np.max(np.abs((u*np.conj(v)-v*np.conj(u))/1j-1)))
            ok=err<=1e-9 and (ex<=1e-9 or rel<=1e-4) and wr<=1e-8
            cases.append(dict(k=k,eps=eps,occupation_numeric=num,occupation_exact=ex,
                              absolute_error=err,relative_error=rel,wronskian_max=wr,
                              nfev=sol.nfev,passed=ok))
    ks=np.array([8.,16.,32.,64.,128.,256.])
    wm=np.sqrt(ks*ks+1);wp=np.sqrt(ks*ks+4)
    sudden=9/(4*wm*wp*(wp+wm)**2);rat=ks**4*sudden/(9/16)
    grid=np.geomspace(.5,256,500)
    wm=np.sqrt(grid*grid+1);wp=np.sqrt(grid*grid+4)
    curve=dict(k=grid.tolist(),sudden=(9/(4*wm*wp*(wp+wm)**2)).tolist(),
               smooth={str(e):spectrum(grid,e).tolist() for e in [.1,.3,1.]})
    passed=all(c["passed"] for c in cases) and abs(rat[-1]-1)<1e-3
    out=dict(status="PASS" if passed else "FAIL",cases=cases,
       abrupt=dict(k=ks.tolist(),occupation=sudden.tolist(),scaled_k4_ratio=rat.tolist(),
                   final_relative_error=float(abs(rat[-1]-1)),asymptote=9/16),
       source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       scope="Exact and numerical toy controls; not archived high-k evolution")
    Path("outputs").mkdir(exist_ok=True)
    Path("outputs/smooth_controls.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    Path("outputs/smooth_spectra.json").write_text(json.dumps(curve,separators=(",",":"),allow_nan=False)+"\n")
    print("HDBLAST_FRW_SMOOTH_JSON="+json.dumps(out,separators=(",",":"),allow_nan=False))
    return 0 if passed else 1
if __name__=="__main__":raise SystemExit(run())

