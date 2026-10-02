"""Post hoc smooth-width limit of the exact toy occupation, not archived energy."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from smooth_controls import spectrum
def energy(eps,refine=False):
    # log-k integrand including dk=k dlogk. Mass in the final static region is 2.
    lo=-18. if refine else -15.
    hi=np.log((70. if refine else 50.)/eps)
    def f(y):
        k=np.exp(y)
        return float(k**3*np.sqrt(k*k+4)*spectrum(k,eps)/(2*np.pi**2))
    val,err=quad(f,lo,hi,epsabs=2e-11 if refine else 1e-9,epsrel=2e-11 if refine else 1e-9,limit=300)
    return float(val),float(err)
cases=[]
for e in [.1,.03,.01,.003,.001,.0003]:
    E,err=energy(e);Et,errt=energy(e,True)
    cases.append(dict(epsilon=e,energy=E,quadrature_estimate=err,extended_band_tight_energy=Et,
                      refinement_absolute_difference=abs(E-Et)))
C=9/(32*np.pi**2);slopes=[]
for a,b in zip(cases[:-1],cases[1:]):
    slope=(b["energy"]-a["energy"])/np.log(a["epsilon"]/b["epsilon"])
    slopes.append(dict(from_epsilon=a["epsilon"],to_epsilon=b["epsilon"],logarithmic_slope=float(slope),
                       leading_coefficient=float(C),relative_difference=float((slope-C)/C)))
report=dict(status="EXECUTED_POST_HOC_EXACT_SPECTRUM_QUADRATURE",cases=cases,slopes=slopes,
    analytic_prediction="E(epsilon)=DeltaU^2/(32*pi^2)*ln(1/epsilon)+O(1) as epsilon tends to zero, toy DeltaU=3,a=1.",
    scope="Late particle energy in the separate static-ended toy step. No HDBLAST energy or fitted physical smoothing width.",
    caution="Quadrature estimates/refinement differences are not rigorous interval enclosures; no registered threshold changed.")
Path("outputs/smoothing_limit.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
print("HDBLAST_FRW_SMOOTHING_LIMIT_JSON="+json.dumps(report,separators=(",",":"),allow_nan=False))

