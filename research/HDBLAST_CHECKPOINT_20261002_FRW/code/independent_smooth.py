"""Independent four-real-component Radau check; no production-code imports."""
import json, math, hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
def logsinh(x):
    if x<=0: raise ValueError("positive argument required")
    return math.log(math.sinh(x)) if x<20 else x-math.log(2)+math.log1p(-math.exp(-2*x))
def exact(k,eps):
    wm=math.sqrt(k*k+1);wp=math.sqrt(k*k+4)
    return math.exp(2*logsinh(math.pi*eps*3/(wp+wm)/2)-logsinh(math.pi*eps*wp)-logsinh(math.pi*eps*wm))
def case(k,eps):
    wm=math.sqrt(k*k+1);wp=math.sqrt(k*k+4)
    ti,tf=-18*eps,18*eps
    step=min(eps/20,.2/wp)
    u0=np.exp(-1j*wm*ti)/math.sqrt(2*wm);v0=-1j*wm*u0
    def freq(t):return k*k+1+1.5*(1+math.tanh(t/eps))
    def rhs(t,y):
        w2=freq(t);return np.array([y[2],y[3],-w2*y[0],-w2*y[1]])
    def jac(t,y):
        w2=freq(t);return np.array([[0,0,1,0],[0,0,0,1],[-w2,0,0,0],[0,-w2,0,0]])
    sol=solve_ivp(rhs,(ti,tf),[u0.real,u0.imag,v0.real,v0.imag],
        method="Radau",jac=jac,rtol=2e-11,atol=2e-13,max_step=step)
    if not sol.success:raise RuntimeError(sol.message)
    u=sol.y[0]+1j*sol.y[1];v=sol.y[2]+1j*sol.y[3]
    beta=np.exp(-1j*wp*tf)*(math.sqrt(wp/2)*u[-1]-1j*v[-1]/math.sqrt(2*wp))
    alpha=np.exp(1j*wp*tf)*(math.sqrt(wp/2)*u[-1]+1j*v[-1]/math.sqrt(2*wp))
    n=float(abs(beta)**2);ne=exact(k,eps);ae=abs(n-ne);re=ae/ne
    wr=float(np.max(np.abs(u*np.conj(v)-v*np.conj(u)-1j)))
    checks=dict(absolute=ae<=1e-9,relative_if_applicable=ne<=1e-9 or re<=1e-4,wronskian=wr<=1e-8)
    return dict(k=k,eps=eps,occupation_numeric=n,occupation_exact=ne,absolute_error=ae,
       relative_error=re,relative_gate_applies=ne>1e-9,wronskian_max=wr,
       bogoliubov_error=float(abs(abs(alpha)**2-abs(beta)**2-1)),
       tail_L1_each=1.5*eps*math.log1p(math.exp(-36)),
       tail_note="frequency-squared tail integral, not a certified state-error bound",
       nfev=sol.nfev,njev=sol.njev,nlu=sol.nlu,nodes=len(sol.t),checks=checks,passed=all(checks.values()))
def main():
    cases=[case(k,e) for e in (.1,.3,1.) for k in (.5,2.,8.)]
    out=dict(status="PASS" if all(c["passed"] for c in cases) else "FAIL",
       solver="Radau, four real components, analytic Jacobian",rtol=2e-11,atol=2e-13,
       source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=cases,
       scope="Smooth oscillator control only; no archived renormalized source")
    Path("outputs").mkdir(exist_ok=True)
    Path("outputs/independent_smooth.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("HDBLAST_FRW_RADAU_JSON="+json.dumps(out,separators=(",",":"),allow_nan=False))
    return 0 if out["status"]=="PASS" else 1
if __name__=="__main__":raise SystemExit(main())

