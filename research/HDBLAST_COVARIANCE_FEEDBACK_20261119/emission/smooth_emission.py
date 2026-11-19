"""Post-plan quiet readout for the two preparations matching initial derivatives."""
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


def summarize(values):
    summary={}
    for win in pairs.WINDOWS:
        for band in ("full","screen"):
            row={}
            def sample(wall,prep,key):
                return values[wall][prep][win][band][observable][key]
            for observable in ("number","energy"):
                for key in ("total","excess","total_density"):
                    ref=sample("erase","erase","total_density" if key=="total_density" else "total")
                    fixed=sample("erase","smoothplus",key)-sample("erase","smoothminus",key)
                    responding=sample("smoothplus","smoothplus",key)-sample("smoothminus","smoothminus",key)
                    midpoint=(sample("smoothplus","smoothplus",key)+sample("smoothminus","smoothminus",key)
                              -2*sample("erase","erase",key))
                    linear=(sample("erase","smoothplus",key)+sample("erase","smoothminus",key)
                            -2*sample("erase","erase",key))
                    row[observable+"_"+key]=dict(
                        state_difference_fixed=fixed/ref, state_difference_responding=responding/ref,
                        geometry_difference=(responding-fixed)/ref,
                        midpoint_fixed=linear/ref,midpoint_responding=midpoint/ref)
            summary.setdefault(win,{})[band]=row
    return summary


def run(repo,native,supplement):
    method=pairs.published_method(repo)
    sfork=json.loads((supplement/"FORK.json").read_bytes())
    mfork=json.loads((native/"FORK.json").read_bytes())
    out={}
    with zipfile.ZipFile(native/"WALLS.zip") as previous, zipfile.ZipFile(supplement/"WALLS.zip") as current:
        for grid,fork in sfork.items():
            q,weights,states=main.fork_states(fork)
            qm,wm,refstates=main.fork_states(mfork[grid])
            if not np.array_equal(q,qm) or not np.array_equal(weights,wm):
                raise ValueError("Native momentum rule changed")
            states["erase"]=refstates["erase"]
            primary,other={},{}
            for wall in ("retain","erase","smoothplus","smoothminus"):
                archive=previous if wall in ("retain","erase") else current
                with archive.open(f"{grid}/{wall}/quiet/wall.csv") as stream:
                    d=main.csv_data(stream)
                start,bg=pairs.wall_interpolants(d,fork["reference_A"])
                af=float(np.exp(bg[1](pairs.CUT)))
                coeff=pairs.fourier_response(start,bg,q,1024,32)
                check=pairs.direct_response(method,bg,q,16384,16)
                primary[wall]={p:pairs.evaluate(coeff,n,c,weights,af) for p,(n,c) in states.items()}
                other[wall]={p:pairs.evaluate(check,n,c,weights,af) for p,(n,c) in states.items()}
            first,second=summarize(primary),summarize(other)
            max_difference=max(abs(value-second[win][band][obs][key]) for win,bands in first.items()
                for band,observables in bands.items() for obs,values in observables.items()
                for key,value in values.items())
            out[grid]=dict(values=primary,normalized_differences=first,
                           independent_differences=second,max_difference_between_methods=max_difference)
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    for key in ("repository","native","native_smooth"):ap.add_argument(key,type=Path)
    args=ap.parse_args()
    print(json.dumps(run(args.repository.resolve(),args.native,args.native_smooth),indent=2,allow_nan=False))
