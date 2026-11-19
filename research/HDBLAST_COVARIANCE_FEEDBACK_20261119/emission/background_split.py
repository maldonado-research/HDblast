"""Separate the scale-factor and scalar-mass contributions to a quiet readout."""
import argparse
import importlib.util
import json
from pathlib import Path
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("main_emission",HERE/"emission.py")
main=importlib.util.module_from_spec(spec)
spec.loader.exec_module(main)
pairs=main.pairs


def split(values):
    ref=values["retain/retain"]
    a=values["reflect/retain"]
    m=values["retain/reflect"]
    both=values["reflect/reflect"]
    out={}
    for obs,row in ref.items():
        out[obs]={}
        for key in ("total","total_density","excess"):
            da=.5*((a[obs][key]-row[key])+(both[obs][key]-m[obs][key]))
            dm=.5*((m[obs][key]-row[key])+(both[obs][key]-a[obs][key]))
            cross=both[obs][key]-a[obs][key]-m[obs][key]+row[key]
            scale=row["total_density" if key=="total_density" else "total"]
            out[obs][key]={name:float(value/scale) for name,value in (
                ("scale_factor",da),("mass",dm),("interaction",cross),
                ("total",both[obs][key]-row[key]))}
    return out


def run(repo,native):
    method=pairs.published_method(repo)
    forks=json.loads((native/"FORK.json").read_bytes())
    out={}
    with zipfile.ZipFile(native/"WALLS.zip") as walls:
        for grid,f in forks.items():
            q,weight,states=main.fork_states(f)
            n,c=states["reflect"]
            bg={}
            starts={}
            for prep in ("retain","reflect"):
                with walls.open(f"{grid}/{prep}/quiet/wall.csv") as stream:
                    data=main.csv_data(stream)
                starts[prep],bg[prep]=pairs.wall_interpolants(data,f["reference_A"])
            if abs(starts["retain"]-starts["reflect"]) > 2e-12:
                raise ValueError("Fork mismatch")
            primary={}
            independent={}
            for a, (clock,loga,_,elapsed) in bg.items():
                for m, (_,_,mass,_) in bg.items():
                    history=(clock,loga,mass,elapsed)
                    af=float(np.exp(loga(pairs.CUT)))
                    first=pairs.fourier_response(starts[a],history,q,1024,32)
                    other=pairs.direct_response(method,history,q,16384,16)
                    primary[a+"/"+m]=pairs.evaluate(first,n,c,weight,af)["wide"]["full"]
                    independent[a+"/"+m]=pairs.evaluate(other,n,c,weight,af)["wide"]["full"]
            x,y=split(primary),split(independent)
            uncertainty=max(abs(value-y[obs][k][field]) for obs,rows in x.items()
                            for k,v in rows.items() for field,value in v.items())
            out[grid]={"four_corners":primary,"symmetric_contributions":x,
                       "independent_contributions":y,"max_method_difference":uncertainty}
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("repository","native"):ap.add_argument(key,type=Path)
    args=ap.parse_args()
    print(json.dumps(run(args.repository.resolve(),args.native),indent=2,allow_nan=False))
