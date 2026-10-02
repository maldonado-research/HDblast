"""Registered selected-mode primary/independent comparison; requires protocol hash.

Only the comparison harness imports the primary DOP853 implementation. The
independent Radau implementation has no production-code imports.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


def import_file(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--protocol",type=Path,required=True)
    parser.add_argument("--protocol-sha256",required=True)
    parser.add_argument("--primary",type=Path,required=True)
    parser.add_argument("--independent",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    actual=hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    if actual!=args.protocol_sha256:
        raise ValueError("Protocol hash mismatch; no field evolution started")
    independent=json.loads(args.independent.read_text())
    if independent["protocol_sha256"]!=actual:
        raise ValueError("Independent data were generated under a different protocol")
    ownpath=Path(__file__).with_name("independent_radau.py")
    own=import_file(ownpath,"independent_mode_reference")
    if independent["source_sha256"]!=hashlib.sha256(ownpath.read_bytes()).hexdigest():
        raise ValueError("Independent implementation differs from data provenance")
    primary=import_file(args.primary,"primary_selected_mode_comparison")
    times=np.asarray(own.TIMES)
    momenta=np.asarray(own.MOMENTA)
    cases=[]
    for amplitude in own.AMPLITUDES:
        u,v,work,solver=primary.evolve(momenta,times,amplitude,primary.TIGHT)
        for index,k in enumerate(momenta):
            reference=next(c for c in independent["cases"] if c["A"]==amplitude and c["k"]==k)
            samples=[]
            for ti,eta in enumerate(times):
                expected=reference["samples"][ti]
                if expected["eta"]!=eta:
                    raise ValueError("Time node mismatch")
                ue=complex(*expected["u"])
                ve=complex(*expected["v"])
                obs=own.observables(float(eta),u[ti,index],v[ti,index],float(k),amplitude)
                gates={"u":dict(difference=abs(u[ti,index]-ue),limit=2e-8),
                       "v":dict(difference=abs(v[ti,index]-ve),limit=2e-8)}
                for key in ("rho","p","Q"):
                    scale=max(abs(obs[key]),abs(expected[key]))
                    gates[key]=dict(difference=abs(obs[key]-expected[key]),
                                    limit=2e-7+2e-8*scale)
                for gate in gates.values():
                    gate["difference"]=float(gate["difference"])
                    gate["passed"]=bool(gate["difference"]<=gate["limit"])
                samples.append(dict(eta=float(eta),
                    primary_u=[float(u[ti,index].real),float(u[ti,index].imag)],
                    primary_v=[float(v[ti,index].real),float(v[ti,index].imag)],
                    primary_bare={key:float(obs[key]) for key in ("rho","p","Q")},gates=gates))
            cases.append(dict(A=amplitude,k=float(k),solver=solver,samples=samples,
                passed=all(g["passed"] for row in samples for g in row["gates"].values())))
    result=dict(status="PASS" if all(c["passed"] for c in cases) else "FAIL",
        scope="Registered selected-mode independent solver comparison; no momentum integral",
        protocol_sha256=actual,
        harness_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        primary_source_sha256=hashlib.sha256(args.primary.read_bytes()).hexdigest(),
        independent_source_sha256=independent["source_sha256"],
        independent_result_sha256=hashlib.sha256(args.independent.read_bytes()).hexdigest(),
        future_phase_note="Radau propagates the static future exactly with phase origin eta=T; DOP853 evolves through it. The opposite phases in positive/negative-frequency components affect complex alpha,beta but not occupation. Raw u,v share incoming phase origin eta=0 and are compared directly.",
        cases=cases)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    maximum={key:max(row["gates"][key]["difference"] for c in cases for row in c["samples"])
             for key in ("u","v","rho","p","Q")}
    print(json.dumps(dict(status=result["status"],maximum_absolute_differences=maximum)))
    return 0 if result["status"]=="PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())
