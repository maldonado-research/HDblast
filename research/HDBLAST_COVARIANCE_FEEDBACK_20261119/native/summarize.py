"""Package and compare self-consistently evolved histories on a common proper-time cut."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np

GRIDS=("coarse","fine")
STATES=("retain","erase","reflect","kick")
ARMS=("quiet","plus","minus","far")
CUT=2.10
READOUT=("H_over_H0","phi_b","A_b","v_over_H0","rho_hat","p_hat","j_hat","comoving_number",
         "scalar_work","pressure_work","comoving_energy")


def json_text(obj):
    return json.dumps(obj,indent=2,allow_nan=False)+"\n"


def store(archive,name,content):
    item=zipfile.ZipInfo(name)
    item.compress_type=zipfile.ZIP_DEFLATED
    archive.writestr(item,content)


def package(work,out):
    forks={}
    matches={}
    with zipfile.ZipFile(out/"WALLS.zip","w") as archive:
        for grid in GRIDS:
            forks[grid]=json.loads((work/grid/"FORK.json").read_bytes())
            for state in STATES:
                for arm in ARMS:
                    name=f"{grid}/{state}/{arm}"
                    store(archive,name+"/wall.csv",(work/name/"wall.csv").read_bytes())
                    if arm!="quiet":
                        matches[name]=json.loads((work/name/"match.json").read_bytes())
        store(archive,"matches.json",json_text(matches).encode())
    (out/"FORK.json").write_text(json_text(forks))


def analyze(out):
    forks=json.loads((out/"FORK.json").read_bytes())
    rows={}
    phase_differences={}
    pulse_differences={}
    with zipfile.ZipFile(out/"WALLS.zip") as archive:
        matches=json.loads(archive.read("matches.json"))
        for grid in GRIDS:
            fork=forks[grid]
            for state in STATES:
                for arm in ARMS:
                    name=f"{grid}/{state}/{arm}"
                    d=np.genfromtxt(io.BytesIO(archive.read(name+"/wall.csv")),names=True,delimiter=",")
                    clock=d["H0tau"]
                    if clock[0]>=CUT or clock[-1]<CUT or np.any(np.diff(clock)<=0):
                        raise ValueError("Proper-time coverage "+name)
                    get=lambda key: float(np.interp(CUT,clock,d[key]))
                    cut={key:get(key) for key in READOUT}
                    cross=np.flatnonzero(np.diff(np.signbit(d["phi_b"]-.5)))
                    if not len(cross):
                        raise ValueError("No mass crossing")
                    i=cross[0]
                    u=(.5-d["phi_b"][i])/(d["phi_b"][i+1]-d["phi_b"][i])
                    cut["cross_s"]=float(clock[i]+u*(clock[i+1]-clock[i]))
                    used=clock<=CUT
                    ledger=d["comoving_energy"]-d["scalar_work"]+d["pressure_work"]
                    err=float(np.max(abs(ledger[used]-ledger[0])))
                    cut["max_energy_ledger_residual"]=err
                    cut["max_shell_constraint"]={key:float(np.nanmax(d[key][used])) for key in ("Hshell","Mshell")}
                    cut["max_near_constraint"]={key:float(np.nanmax(d[key][used])) for key in ("Hmax_near","Mmax_near")}
                    cut["min_determinant"]=float(np.min(d["min_uncertainty"][used]))
                    cut["min_occupation"]=float(np.min(d["min_occupation"][used]))
                    if abs(fork["H0tau"]-clock[0])>1.e-10 or cut["min_determinant"]<.24999 or err>1.e-3:
                        raise ValueError("Initial cut, physicality, or energy ledger mismatch "+name)
                    rows[name]=cut
                    if state!="retain":
                        ref=rows[f"{grid}/retain/{arm}"]
                        phase_differences[name]={k:cut[k]-ref[k] for k in (*READOUT,"cross_s")}
                    if arm!="quiet":
                        ref=rows[f"{grid}/{state}/quiet"]
                        pulse_differences[name]={k:cut[k]-ref[k] for k in (*READOUT,"cross_s")}
    pulse_max={k:max(abs(r[k]) for r in matches.values()) for k in
               ("matched_delta_y_left","H_relative_delta_max","M_relative_delta_max")}
    outdata=dict(cut=CUT,rows=rows,phase_minus_retained=phase_differences,pulse_minus_quiet=pulse_differences,
                 shooting_residual_max=pulse_max)
    (out/"NATIVE.json").write_text(json_text(outdata))
    return outdata


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder",type=Path,help="Archive/report destination")
    ap.add_argument("--from-work",type=Path,help="Native scratch directory to package before analysis")
    args=ap.parse_args()
    args.folder.mkdir(parents=True,exist_ok=True)
    if args.from_work:
        package(args.from_work,args.folder)
    data=analyze(args.folder)
    print(f"Compared {len(data['rows'])} native solutions at s={CUT}")


if __name__=="__main__":main()
