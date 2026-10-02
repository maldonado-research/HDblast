"""Postprocess saved runs; no further mode integration."""
import json
import os
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
import pv_crossing as p
base=Path("outputs/pv_actual")
r=json.loads((base/"pv_actual.json").read_text())
cases={x["name"]:x for x in r["cases"]}
curves={name:np.load(base/(name+"_curves.npz"),allow_pickle=False)["curves"] for name in cases}
primary=curves["primary"]
def direct_difference(c):
    return {key:float(np.max(np.abs(c[:,col]-primary[:,col]))/scale)
            for key,col,scale in [("rho",9,p.M**4),("p",10,p.M**4),("S",11,p.M**2)]}
direct={name:direct_difference(c) for name,c in curves.items() if name!="primary"}
h=p.COMMON[1]-p.COMMON[0]
trace={}
finite_trace={}
for name,c in curves.items():
    Q=2*c[:,11]
    local=(c[:,1]-p.R)**2/(32*np.pi**2)
    local+=np.array([p.pulse(t)[2] for t in c[:,0]])/(96*np.pi**2)
    derivative=.5*p.second(Q,h,9)-c[:,1]*Q
    expected=derivative+local
    target=-c[:,9]+3*c[:,10]
    trace[name]=float(np.nanmax(np.abs(expected[8:-8]-target[8:-8]))/p.M**4)
    K=cases[name]["K"];v=K/np.sqrt(K*K+p.R)
    finite_trace[name]=float(np.nanmax(np.abs((derivative+v**3*local-target)[8:-8]))/p.M**4)
cross={}
for name,c in curves.items():
    z=c[1000]
    cross[name]=dict(t=float(z[0]),rho=float(z[2]),p=float(z[3]),S=float(z[4]),
                    rho_direct=float(z[9]),p_direct=float(z[10]),S_direct=float(z[11]),
                    p_potential_only=float(z[12]))
negative=[]
for a,b in [("lambda4","lambda8"),("lambda8","primary")]:
    la=cases[a]["lambda_"];lb=cases[b]["lambda_"]
    prediction=-p.pulse(0.)[2]*np.log(lb/la)/(48*np.pi**2)
    observed=cross[b]["p_potential_only"]-cross[a]["p_potential_only"]
    negative.append(dict(from_lambda=la,to_lambda=lb,observed=observed,
                         asymptotic_prediction=float(prediction),
                         fractional_discrepancy=float((observed-prediction)/abs(prediction))))
z=np.load(base/"primary_curves.npz",allow_pickle=False)
k=z["k"];w=z["weights"];n=z["n_exact"][0]
number=float(w@n);energy=float(w@(np.sqrt(k*k+p.R)*n))
gas_p=float(w@(k*k/(3*np.sqrt(k*k+p.R))*n))
out=dict(status="EXECUTED_POSTPROCESSING_OF_REGISTERED_DATA",
         first_reference_run=36949788139,execution_run=os.environ.get("GITHUB_RUN_ID"),source_ref=r["source_ref"],
         direct_cutoff_and_numeric_differences=direct,
         direct_limit_trace_residual_scaled=trace,
         direct_finite_K_trace_residual_scaled=finite_trace,
         finite_K_trace_note="Post hoc exact factor v_r^3 multiplies the local matched trace term; Q=2S. No registered threshold changed.",
         direct_trace_note="Continuum matched trace used against finite-K direct integrals; residual includes the unremoved cutoff tail.",
         crossing=cross,pressure_negative_control=negative,
         asymptotic_exact_particles=dict(number=number,energy=energy,gas_pressure=gas_p,
                         gas_w=gas_p/energy,mean_energy=energy/number,
                         endpoint_rho=cases["primary"]["endpoint_rho"],
                         endpoint_energy_relative_difference=(cases["primary"]["endpoint_rho"]-energy)/energy),
         caution="Natural-scale error targets are not relative errors on the small source. Finite-Lambda convergence is not proof of a regulator limit.")
Path("outputs/analysis_summary.json").write_text(json.dumps(out,indent=2)+"\n")
print("HDBLAST_PV_ANALYSIS_JSON="+json.dumps(out,separators=(",",":")))
