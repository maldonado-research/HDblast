"""Read-only audit of full replay arrays; no field evolution.

The exact previously executed helper is retained in history/audit_matrix_executed.py.
This public CLI adds a required case-count guard and omits absolute input paths
from the report. All numerical audit expressions and tolerances are unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path,help="Completed full runtime matrix directory")
    parser.add_argument("--expected-cases",type=int,choices=(4,13),required=True,
                        help="13 for the original matrix; 4 for the separate follow-up")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    summary=json.loads((args.directory/"summary.json").read_text())
    if summary["status"] not in ("PASS","FAIL"):
        raise ValueError("Matrix is incomplete; no final audit produced")
    paths=sorted(args.directory.glob("*/diagnostics.json"))
    if len(paths)!=args.expected_cases:
        raise ValueError("Missing or unexpected completed cases; no final audit produced")
    cases=[]
    for path in paths:
        diagnostics=json.loads(path.read_text())
        with np.load(path.with_name("arrays.npz"),allow_pickle=False) as archive:
            data={key:archive[key] for key in archive.files}
        with np.load(path.with_name("physical_modes.npz"),allow_pickle=False) as archive:
            for key in archive.files:
                if key in data and not np.array_equal(data[key],archive[key]):
                    raise ValueError("Physical/observable archives disagree on "+key)
                data[key]=archive[key]
        finite=all(np.all(np.isfinite(value)) for value in data.values())
        finite=finite and all(value.dtype!=object for value in data.values())
        t,k,weights,u,v,a=[data[key] for key in ("t","k","weights","u","v","a")]
        wr=float(np.max(abs((u*v.conjugate()-v*u.conjugate())/1j-1)))
        record=dict(case=path.parent.name,finite_numeric_arrays=finite,
                    stored_wronskian_matches=abs(wr-diagnostics["wronskian_max"])<=1e-15,
                    wronskian_max_recomputed=wr,cutoffs={})
        for key,reported in diagnostics["cutoffs"].items():
            K=int(key)
            mask=k<K
            errors={}
            for quantity in ("rho","p","Q"):
                total=np.sum(np.asarray(data["mode_"+quantity][:,mask],dtype=np.longdouble)
                             *np.asarray(weights[mask],dtype=np.longdouble)[None,:],axis=1,
                             dtype=np.longdouble)
                errors[quantity]=float(np.max(abs(total-data[f"K{K}_{quantity}"])))
            future=(t>=1) if not diagnostics["static"] else np.ones(len(t),dtype=bool)
            wf=np.sqrt(k*k+a[-1]**2*4.)
            beta=(wf[None,:]*u[future]-1j*v[future])/np.sqrt(2*wf[None,:])
            energy=np.sum((abs(beta[:,mask])**2*wf[mask][None,:]/a[-1]**4)
                          *weights[mask][None,:],axis=1)
            stored_occ=data[f"K{K}_future_occ_rho"]
            occerror=float(np.max(abs(energy-stored_occ)))
            scalar_gate=reported["future_spectral_energy"]
            actualdiff=float(np.max(abs(data[f"K{K}_rho"][future]-energy)))
            record["cutoffs"][key]=dict(reintegration_errors=errors,
                reintegration_passed=max(errors.values())<1e-10,
                future_occupation_recompute_error=occerror,
                future_occupation_recompute_passed=occerror<1e-10,
                future_spectral_absolute_difference=actualdiff,
                reported_spectral_difference_matches=abs(actualdiff-scalar_gate["absolute_difference"])<1e-10)
        record["passed"]=finite and record["stored_wronskian_matches"] and all(
            c["reintegration_passed"] and c["future_occupation_recompute_passed"]
            and c["reported_spectral_difference_matches"] for c in record["cutoffs"].values())
        cases.append(record)
    report=dict(scope="Post-execution independent read-only artifact audit, no new field evolution",
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                run_name=args.directory.name,expected_cases=args.expected_cases,
                matrix_status=summary["status"],
                audit_status="PASS" if all(c["passed"] for c in cases) else "FAIL",cases=cases)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(audit_status=report["audit_status"],matrix_status=summary["status"],cases=len(cases))))
    return 0 if report["audit_status"]=="PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())
