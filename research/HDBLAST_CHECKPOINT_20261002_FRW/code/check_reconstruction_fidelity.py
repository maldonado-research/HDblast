"""Post hoc diagnostics of interpolation and the native high-precision evaluator.
No new physical acceptance threshold; executed separately from registered verdict.
"""
import json
from pathlib import Path
import numpy as np
import geometry_audit as ga
base=Path(__file__).resolve().parent.parent
prior=ga.read_prior(base/"code/legacy_quantum_source.py")
cases=[]
for tag in (prior.FINE,prior.COARSE):
    data,_=prior.load(base/"data/A1",tag)
    saved=data["H0tau"];mask=(saved>=0)&(saved<=ga.END);tt=saved[mask]
    for kind in ga.KINDS:
        la,ph=ga.build(prior,data,kind)
        nodal=dict(ln_a_max_abs=float(np.max(np.abs(la(tt)-(data["ln_a"][mask]-data["ln_a"][0])))),
                   phi_max_abs=float(np.max(np.abs(ph(tt)-data["phi_b"][mask]))))
        native=[]
        if isinstance(la,ga.NativeSmoothSpline):
            knots=np.unique(la.x);interior=knots[(knots>0)&(knots<ga.END)]
            probes=np.unique(np.concatenate([np.linspace(0.,ga.END,33),interior[::17]]))
            for label,poly in [("ln_a",la),("phi",ph)]:
                for order in range(poly.bs.k+1):
                    bs=np.array([float(poly(t,nu=order)) for t in probes])
                    hp=np.array([float(poly.one_sided_mp(float(t),False,order)) for t in probes])
                    if not (np.isfinite(bs).all() and np.isfinite(hp).all()):
                        raise ValueError("Nonfinite evaluator comparison")
                    err=float(np.max(np.abs(bs-hp)))
                    native.append(dict(field=label,order=order,n_probes=len(probes),max_abs=err,
                                       scaled_by_sampled_max=err/max(1.,float(np.max(np.abs(hp))))))
        cases.append(dict(tag=tag,kind=kind,nodal_value_fidelity=nodal,native_vs_80digit_deBoor=native))
report=dict(status="EXECUTED_POST_HOC_REPRESENTATION_DIAGNOSTIC",cases=cases,
    scope="Stored-value interpolation and two numerical evaluations of the same spline; not physical derivative accuracy or a new registered criterion.")
out=base/"outputs"/"reconstruction_fidelity.json"
out.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
print("HDBLAST_FRW_FIDELITY_JSON="+json.dumps(report,separators=(",",":"),allow_nan=False))
