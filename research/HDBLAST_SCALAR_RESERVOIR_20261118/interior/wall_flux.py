"""Recompute the wall boundary comparison from both derived B1 CSV files."""
import csv
import hashlib
import io
import json
import math
from bisect import bisect_left
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def integral(times, values, absolute=False):
    pieces=[]
    for a,b,l,r in zip(times,times[1:],values,values[1:]):
        if absolute and l*r<0:
            pieces.append((b-a)*(l*l+r*r)/(2*(abs(l)+abs(r))))
        else:
            term=(b-a)*(l+r)/2
            pieces.append(abs(term) if absolute else term)
    return math.fsum(pieces)


def value(t,v,x):
    i=bisect_left(t,x)
    if t[i]==x:
        return v[i]
    f=(x-t[i-1])/(t[i]-t[i-1])
    return v[i-1]+f*(v[i]-v[i-1])


def bin_summary(t,rate,lo,hi):
    x=[lo]+[v for v in t if lo<v<hi]+[hi]
    return integral(x,[value(t,rate,p) for p in x])


def analyze(grid, receipt):
    data=(ROOT/f"wall_flux_{grid}.csv").read_bytes()
    if hashlib.sha256(data).hexdigest()!=receipt["csv_sha256"]:
        raise ValueError("CSV hash")
    rows=list(csv.DictReader(io.StringIO(data.decode("ascii"))))
    vals={k:[float(r[k]) for r in rows] for k in rows[0]}
    t=vals["s"]
    if len(t)!=155 or t[0]!=1 or t[-1]!=4 or not all(y>x for x,y in zip(t,t[1:])):
        raise ValueError("wrong time samples")
    bc=vals["junction_rate"]
    sig=integral(t,vals["tension_rate"])
    mat=integral(t,vals["matter_rate"])
    law=integral(t,bc)
    if abs(sig+mat-law)>1e-10:
        raise ValueError("source decomposition")
    out=dict(csv_sha256=receipt["csv_sha256"], tension=sig,matter=mat,junction=law,
             methods={})
    for order in range(1,5):
        for label in ("interior","with_wall"):
            k=f"{label}_{order}"
            measured=vals[k]
            area=integral(t,measured)
            err=[a-b for a,b in zip(measured,bc)]
            out["methods"][k]=dict(integral=area,inferred_matter_after_tension=area-sig,
                signed_flux_minus_junction=integral(t,err),signed_difference_over_matter=integral(t,err)/mat,
                absolute_difference_over_absolute_junction=integral(t,err,True)/integral(t,bc,True),
                absolute_difference_over_absolute_matter=integral(t,err,True)/integral(t,vals["matter_rate"],True))
    key=vals["interior_4"]
    out["fourth_order_interior_bins"]=[]
    for i in range(8):
        lo,hi=1+3*i/8,1+3*(i+1)/8
        flux=bin_summary(t,key,lo,hi)
        target=bin_summary(t,bc,lo,hi)
        out["fourth_order_interior_bins"].append(dict(interval=[lo,hi],flux=flux,junction=target,
                flux_minus_junction=flux-target, matter=bin_summary(t,vals["matter_rate"],lo,hi)))
    if not law<0 or not all(v["integral"]<0 for v in out["methods"].values()):
        raise ValueError("integrated flux sign")
    return out


if __name__=="__main__":
    rr=json.loads((ROOT/"WALL_RERUN.json").read_text())
    result={grid:analyze(grid,rr[grid]) for grid in ("coarse","fine")}
    print(json.dumps(result,indent=2))
