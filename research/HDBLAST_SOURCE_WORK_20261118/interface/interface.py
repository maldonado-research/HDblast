"""Integrate two-sided bulk scalar flux, shell-tension work and reciprocal matter work."""
import argparse
import bisect
import csv
import io
import json
import math
from pathlib import Path
import zipfile


HERE=Path(__file__).resolve().parent
CUT=2.10
ARMS=("quiet","plus","minus","far")


def lagrange_derivative(z):
    """Weights for the derivative at zero from interior samples of a polynomial."""
    return [math.prod(-xi/(xj-xi) for i,xi in enumerate(z) if i!=j)
            * math.fsum(-1/xi for i,xi in enumerate(z) if i!=j) for j,xj in enumerate(z)]


def interpolate(x,y,t):
    if not x[0] <= t <= x[-1]:raise ValueError("Uncovered proper or conformal clock")
    i=min(len(x)-2,bisect.bisect_right(x,t)-1)
    return y[i]+(y[i+1]-y[i])*(t-x[i])/(x[i+1]-x[i])


def quad3(x,y,lo,hi):
    # Integrate the quadratic Lagrange interpolant in a local centered coordinate.
    p=[(t-lo)/(hi-lo) for t in x]
    terms=[]
    for j in range(3):
        a,b=p[(j+1)%3],p[(j+2)%3]
        weight=(1/3-(a+b)/2+a*b)/((p[j]-a)*(p[j]-b))
        terms.append(weight*y[j])
    return (hi-lo)*math.fsum(terms)


def integral(x,y,t):
    i=bisect.bisect_right(x,t)-1
    even=i-i%2
    parts=[quad3(x[j:j+3],y[j:j+3],x[j],x[j+2]) for j in range(0,even,2)]
    if t>x[even]:
        j=min(len(x)-3,even)
        parts.append(quad3(x[j:j+3],y[j:j+3],x[even],t))
    return math.fsum(parts)


def trapezoid(x,y,t):
    i=bisect.bisect_right(x,t)-1
    total=math.fsum((x[j+1]-x[j])*(y[j]+y[j+1])/2 for j in range(i))
    if t>x[i]:
        total+=(t-x[i])*(y[i]+interpolate(x,y,t))/2
    return total


def measure(samples,params,grid):
    d={k:[float(row[k]) for row in samples] for k in samples[0]}
    t,s=d['T'],d['H0tau']
    cut=interpolate(s,t,CUT)
    rb,b=params["rb"][grid],params["b"]
    v=[math.exp(3*(a-params["A_ref"][grid])) for a in d["A_b"]]
    z=params["interior_z"][grid]
    if not all(x<0 for x in z) or any(x!=0 for x in d["gradient_wall_z"]):
        raise ValueError("Expected a fixed boundary and strictly interior stencil")
    phiT=d["bulk_phi_T"]
    for idx in range(len(t)):
        if abs(phiT[idx]-d["v_over_H0"][idx]*math.exp(d["B_b"][idx])/rb)>1e-12:
            raise ValueError("Coordinate velocity mismatch")
        if abs(d["j_hat"][idx]*b/rb-d["J0"][idx])>1e-13:
            raise ValueError("Misnormalized source")
    rates={"tension":[V*pt*sig for V,pt,sig in zip(v,phiT,d["sigma1"])],
           "matter":[V*pt*b/rb*j for V,pt,j in zip(v,phiT,d["j_hat"])]}
    n=len(t)
    for degree in (2,4):
        weights=lagrange_derivative(z[:degree+1])
        grad=[math.fsum(weights[i]*(d[f"inside_phi_{i+1}"][j]-d["inside_phi_1"][j])
                       for i in range(degree+1)) for j in range(n)]
        rates[f"bulk_{degree}"]=[-2*V*math.exp(-B)*pt*pz for V,B,pt,pz in zip(v,d["B_b"],phiT,grad)]
        rates[f"matter_from_bulk_{degree}"]=[a-b for a,b in zip(rates[f"bulk_{degree}"],rates["tension"])]
    w=rates["matter"]
    full={k:integral(t,r,cut) for k,r in rates.items()}
    # Quiet-subtracted excesses are formed only after each history is integrated to the same proper cut.
    rk=b/rb*(interpolate(s,d["scalar_work"],CUT)-d["scalar_work"][0])
    # Independently integrate the shell matter balance in extensive form.
    # The comoving energy and pressure-work states are accumulated separately in the coupled solver.
    energy=b/rb*(interpolate(s,d["comoving_energy"],CUT)-d["comoving_energy"][0])
    pw=b/rb*(interpolate(s,d["pressure_work"],CUT)-d["pressure_work"][0])
    full.update(rk_matter=rk,energy_change=energy,pressure_work=pw,total_matter_signed_defect=energy+pw-rk,
                cut_T=cut,start_T=t[0],start_s=s[0],proper_cut=CUT)
    for degree in (2,4):
        vbulk=full[f"matter_from_bulk_{degree}"]
        full[f"bulk_{degree}_minus_rk_matter"]=vbulk-rk
        full[f"bulk_{degree}_fractional_matter_error"]=(vbulk-rk)/rk
    for step in (1,2,4):
        for quad,rule in (("simpson",integral),("trapezoid",trapezoid)):
            full[f"quartic_{quad}_stride_{step}"]=rule(t[::step],rates["matter_from_bulk_4"][::step],cut)
    full["bulk_rate_max_error_after_tension"]=max(abs(a-b) for tm,a,b in
            zip(t,rates["matter_from_bulk_4"],rates["matter"]) if tm<=cut)
    full["samples"]=len(t)
    full["dt"]=t[1]-t[0]
    full["work_every_other_row"]=integral(t[::2],w[::2],cut)
    full["scalar_work_integral_minus_rk"]=full["matter"]-rk
    full["source_fraction_of_abs_tension"]=rk/abs(full["tension"])
    return full


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samples",type=Path,default=HERE/"SAMPLES.zip",nargs="?")
    parser.add_argument("--walls",type=Path,default=HERE.parent/"coupled/WALLS.zip",
                        help="Baseline coupled record for the cadence cross-check")
    args=parser.parse_args()
    path=args.samples
    result={}
    with zipfile.ZipFile(path) as z, zipfile.ZipFile(args.walls) as walls:
        metadata=json.loads(z.read("parameters.json"))
        for grid in ("coarse","fine"):
            arms={}
            for arm in ARMS:
                with z.open(f"{grid}/{arm}.csv") as f:
                    rows=list(csv.DictReader(io.TextIOWrapper(f,encoding="ascii")))
                stride=2 if grid=="coarse" else 4
                with walls.open(f"{grid}/y1_b1/{arm}.csv") as f:
                    reference=list(csv.DictReader(io.TextIOWrapper(f,encoding="ascii")))
                if len(reference)!=len(rows[::stride]):raise ValueError("Wrong recording cadence")
                errors={k:max(abs(float(a[k])-float(b[k])) for a,b in zip(rows[::stride],reference,strict=True))
                        for k in ("T","H0tau","A_b","phi_b","H_over_H0","j_hat","scalar_work","rho_hat","rad_hat")}
                tolerance={"T":1e-10,"H0tau":1e-10,"A_b":1e-9,"phi_b":1e-9,"H_over_H0":5e-8,
                           "j_hat":5e-4,"scalar_work":1e-4,"rho_hat":1e-4,"rad_hat":1e-4}
                if any(e>tolerance[k] for k,e in errors.items()):raise ValueError("The flux replay changed the wall")
                arms[arm]=measure(rows,metadata,grid)
                arms[arm]["checked_trajectory_errors"]=errors
            quiet=arms["quiet"]
            contrasts={}
            for name in ARMS[1:]:
                row=arms[name]
                sel=("tension","bulk_4","rk_matter","matter","matter_from_bulk_2","matter_from_bulk_4",
                     "energy_change","pressure_work",*(f"quartic_{rule}_stride_{stride}" for rule in
                     ("simpson","trapezoid") for stride in (1,2,4)))
                dx={key:row[key]-quiet[key] for key in sel}
                dx["quartic_error_vs_rk"]=dx["matter_from_bulk_4"]-dx["rk_matter"]
                dx["fractional_error_vs_rk"]=dx["quartic_error_vs_rk"]/dx["rk_matter"]
                dx["fractional_matter_work_change"]=dx["rk_matter"]/quiet["rk_matter"]
                if any(dx[k]*dx["rk_matter"]<=0 for k in ("matter_from_bulk_2","matter_from_bulk_4")):
                    raise ValueError("Incompatible pulse signs")
                limit=.022 if grid=="coarse" else .003
                if abs(dx["fractional_error_vs_rk"])>limit:
                    raise ValueError("Integrated interior-flux reconstruction")
                contrasts[name]=dx
            result[grid]={"arms":arms,"matched_differences":contrasts}
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
