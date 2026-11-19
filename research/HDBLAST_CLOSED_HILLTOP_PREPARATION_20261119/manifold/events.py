"""Resolve radiation equality or turnaround after four finite-velocity extreme launches."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


def load(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def outcome(engine,pot,x,start):
    # Reuse the published two-cycle gate and first crossing 24 units after
    # first minimum, followed by the prior pressureless decay calculation.
    init=engine.constrained(pot,start,x,x+18,0.,2e-10)
    first=init.t_events[0][0]
    switch=engine.gate_time(pot,init)
    stage=engine.constrained(pot,init.sol(switch),switch,first+27,pot.gamma,2e-10)
    cross=next(t for t in stage.t_events[0] if t>=first+24)
    state=engine.values(pot,stage.sol(cross))
    a=state["a"]
    s=state["scalar"]
    out=engine.ENVELOPE.evolve(pot.gamma*a**3*s,1/(a*a*s),state["radiation"]/s)
    return dict(first_event="equality" if out["equality"] else "turnaround",
                ratio_at_maximum=out["turn"]["radiation_over_matter"],
                a_at_first_min=engine.values(pot,init.sol(first))["a"],
                x_first_min=first,x_gate=switch)


def main():
    here=Path(__file__).resolve().parent
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prior",type=Path,default=here.parent/"prior")
    cli.add_argument("--manifold",type=Path,default=here/"MANIFOLD.json")
    cli.add_argument("--output",required=True,type=Path)
    args=cli.parse_args()
    M=load(here/"manifold.py","manifold_solver")
    X=load(here/"extract.py","archive_extractor")
    engine=load(args.prior/"hilltop"/"hilltop.py","hilltop_solver")
    study=json.loads(args.manifold.read_text())
    local=X.extract(args.prior)
    if local["prior_sha256"]!=study["prior_sha256"]:
        raise RuntimeError("Prior input sources different from threshold calculation")
    mod,records=engine.BG.source(args.prior/"archive")
    tables=mod.build_tables()
    result={}
    for k,branch in study["branches"].items():
        field=M.Field(local["branches"][k])
        p=engine.BG.Potential(mod,tables,field.pm,records[k])
        cases=[]
        for i in [0,-1]:
            row=branch["samples"][i]
            base=M.stable(field,row["v"],end=5.5,tol=1e-10)
            for scale in [.94,1.06]:
                seed=p.sign*row["conditional_width"]*scale
                xmatch,point=M.perturbed(field,base,seed,field.row["critical_x"]+.2,branch["match_loga"])
                start=np.array([field.top+point[0],point[1],point[3],0.])
                out=outcome(engine,p,xmatch,start)
                if out["first_event"]!=("equality" if scale<1 else "turnaround"):
                    raise ArithmeticError("Opposite event ordering to the matched bracket")
                cases.append(dict(v=row["v"],offset_multiplier=scale,signed_offset=seed,**out))
                print(k,row["v"],scale,out["first_event"],flush=True)
        result[k]=cases
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
