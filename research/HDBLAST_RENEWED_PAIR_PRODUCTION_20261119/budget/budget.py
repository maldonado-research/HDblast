"""Hamiltonian energy and work for two scalar-clock Yukawa pulses on measured walls."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicHermiteSpline, CubicSpline, PchipInterpolator
from scipy.optimize import brentq

PRIOR = Path("research/HDBLAST_QUANTUM_MEMORY_20261119/source")
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
ORIGINAL = (1.35, 2.05)
LATE = (2.055, 2.145)
READOUT = 2.16
CASES = (("fine","quiet"),("fine","minus"),("coarse","quiet"),("fine","plus"),("fine","far"))
MODELS = {"first":(1.,0.),"second":(0.,1.),"both":(1.,1.),"opposite":(1.,-1.)}
PRIMARY_MODELS = ("first","second","both")
COMPONENTS = ("parent","daughter","interaction","counterterm")

def source_engine(repo, cutoff=120):
    entry=repo/PRIOR/"source.py"
    spec=importlib.util.spec_from_file_location("hdblast_repump_prior_matter",entry)
    engine=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    if cutoff not in (120,240,480):
        raise ValueError("Expected a registered frequency cutoff")
    if cutoff > 120:
        original=engine.make_bath
        def extended(q,high_order,low_order):
            # Add explicit positive oscillators on independent fixed frequency bands.
            om0,g0,_=original(q,high_order,low_order)
            om=[om0]
            g2=[g0]
            x,w=np.polynomial.legendre.leggauss(64)
            left=120
            while left<cutoff:
                right=min(cutoff,2*left)
                omega=left+(right-left)*(1+x)/2
                d_omega=(right-left)*w/2
                om.append(np.broadcast_to(omega,(len(q),len(x))))
                g2.append(omega*(omega**2-q[:,None]**2)*d_omega/(4*np.pi**2))
                left=right
            return np.concatenate(om,axis=1),np.concatenate(g2,axis=1),"piecewise_GL64_to_"+str(cutoff)
        engine.make_bath=extended
        engine.HIGH=cutoff
    return engine

class Gate:
    def __init__(self, quiet_mass, window):
        lo, hi=window
        at=np.linspace(lo,hi,8193)
        phi=quiet_mass(at)/100+.5
        if np.any(np.diff(phi)<=0):
            raise ValueError("Scalar must increase throughout the calibration window")
        self.support=(phi[0],phi[-1])
        self.lo=lo
        self.rate=np.pi/(hi-lo)
        self.inv=CubicSpline(phi,at,bc_type=((1,100/quiet_mass(lo,1)),(1,100/quiet_mass(hi,1))),
                             extrapolate=False)

    def __call__(self, phi):
        active=(phi>self.support[0]) & (phi<self.support[1])
        clipped=np.clip(phi,*self.support)
        s=self.inv(clipped)
        d=self.inv(clipped,1)
        angle=self.rate*(s-self.lo)
        sn,cs=np.sin(angle),np.cos(angle)
        first=self.rate*4*sn**3*cs
        second=4*self.rate**2*(3*sn*sn*cs*cs-sn**4)
        return tuple(np.where(active,z,0.) for z in (sn**4, first*d,
                                                       second*d**2+first*self.inv(clipped,2)))

class Clock:
    def __init__(self, quiet, amplitudes):
        self.gates=(Gate(quiet,ORIGINAL),Gate(quiet,LATE))
        self.amplitudes=amplitudes
        self.support=self.gates[0].support+self.gates[1].support

    def __call__(self, phi):
        primary,late=(g(phi) for g in self.gates)
        a,b=self.amplitudes
        return tuple(a*x+b*y for x,y in zip(primary,late))

def gate_times(engine, native, grid, arm):
    mass=engine.read_native(native,grid,arm)[-1][-1]
    quiet=engine.read_native(native,'fine','quiet')[-1][-1]
    clock=Clock(quiet,(1.,1.))
    return [brentq(lambda s: mass(s)/100+.5-p,mass.x[0],mass.x[-1],xtol=1e-14)
            for p in clock.support]

def make_one(task):
    repo, dest, grid, arm, label, steps, cutoff=task
    engine=source_engine(repo,cutoff)
    engine.ScalarClock=lambda quiet:Clock(quiet,MODELS[label])
    engine.integrate(repo/NATIVE,grid,arm,64,steps=steps,output=dest/label)
    return f"{grid}_{arm}_{label}"

def generate(repo, dest, jobs=2, steps=4096, cases=CASES, models=None, cutoff=120):
    # Each quadratic history uses exactly the same RK graph, including all switch times.
    dest.mkdir(parents=True,exist_ok=True)
    engine=source_engine(repo)
    gates={f"{grid}_{arm}":gate_times(engine,repo/NATIVE,grid,arm) for grid,arm in cases}
    old=dest/"GATES.json"
    known=json.loads(old.read_text()) if old.exists() else {}
    old.write_text(json.dumps(dict(known,**gates),indent=2)+"\n")
    tasks=[(repo,dest,grid,arm,label,steps,cutoff) for grid,arm in cases for label in (models or PRIMARY_MODELS)]
    with ProcessPoolExecutor(max_workers=jobs) as pool:
        for case in pool.map(make_one,tasks):print("Generated",case,flush=True)

def load(tables,grid,arm,scheme,label):
    with np.load(tables/label/f"{grid}_{arm}_N64_{scheme}.npz", allow_pickle=False) as f:
        return {k:f[k] for k in f}

def evaluate(s,y,dy,t):
    return float(CubicHermiteSpline(s,y,dy,extrapolate=False)(t))

def sector(table, start, end):
    s,a,h,vel=(table[k] for k in ("s","a","h","v"))
    if not (s[0] <= start < end <= s[-1]):
        raise ValueError("Interval extends beyond the measured table")
    # Most endpoints are native source knots. Include interpolated endpoints exactly.
    inside=(s>start) & (s<end)
    stamps=np.r_[start,s[inside],end]
    def integral(values):
        curve=CubicSpline(s,values,extrapolate=False)
        return float(simpson(np.r_[curve(start),values[inside],curve(end)],x=stamps))
    energy={}
    for comp in COMPONENTS:
        r=table["rho_"+comp]
        dr=table["drho_"+comp]
        value=a**3*r
        slope=a**3*(dr+3*h*r)
        energy[comp]=[evaluate(s,value,slope,t) for t in (start,end)]
    w_mass=integral(a**3*vel*table["j_parent"])
    w_gate=integral(a**3*vel*(table["j_interaction"]+table["j_counterterm"]))
    gate_power=a**3*vel*(table["j_interaction"]+table["j_counterterm"])
    # All extrema of the cumulative-work interpolant occur at endpoints or at
    # roots of its power. PCHIP preserves each power interval's monotonicity.
    pow_spline=PchipInterpolator(s,gate_power,extrapolate=False)
    roots=pow_spline.roots(extrapolate=False)
    extrema=np.unique(np.r_[start,roots[(roots>start)&(roots<end)],end])
    primitive=pow_spline.antiderivative()
    cumulative=primitive(extrema)-primitive(start)
    increments=np.diff(cumulative)
    gate_path=dict(gross_in=float(np.sum(np.maximum(increments,0.))),
                   gross_out=float(-np.sum(np.minimum(increments,0.))),
                   maximum_cumulative=float(np.max(cumulative)),
                   minimum_cumulative=float(np.min(cumulative)),
                   pchip_net=float(cumulative[-1]))
    w_pressure=integral(-3*h*a**3*table["p"])
    change=sum(pair[1]-pair[0] for pair in energy.values())
    ward=table["drho"]+3*h*(table["rho"]+table["p"])-table["j"]*vel
    residual=change-w_mass-w_gate-w_pressure
    max_ward=float(np.max(abs(ward[(s>=start)&(s<=end)])))
    if abs(residual)>1e-5 or max_ward>1e-8:
        raise ValueError("Energy work or local continuity failed on measured wall")
    return dict(proper_energy=dict((k,b) for k,(_,b) in energy.items()),
                change_energy={k:b-a for k,(a,b) in energy.items()},
                total_proper_change=change, scalar_mass_work=w_mass,
                gate_interaction_work=integral(a**3*vel*table["j_interaction"]),
                gate_completion_work=integral(a**3*vel*table["j_counterterm"]),
                gate_work=w_gate, gate_work_path=gate_path, expansion_work=w_pressure,
                energy_residual=residual, max_pointwise_ward=max_ward,
                peak_abs_current=float(np.max(abs(table["j"][(s>=start)&(s<=end)]))))

def physical_difference(combined, first):
    # Background fields and gate metadata are not matter corrections.
    res=combined.copy()
    for k,v in combined.items():
        if k.startswith(("rho","p_","j","drho","dp","dj","rad","drad","energy","daughter")) or k=="p":
            res[k]=v-first[k]
    return res

def summarize(tables, cases=CASES):
    gates=json.loads((tables/"GATES.json").read_text())
    result={}
    for grid,arm in cases:
        name=grid+"_"+arm
        start=gates[name][2]
        data={"four_gate_knots":gates[name]}
        for scheme in ("static","bare"):
            first,second,both=(load(tables,grid,arm,scheme,label) for label in ("first","second","both"))
            np.testing.assert_array_equal(first['s'],both['s'])
            np.testing.assert_array_equal(first['s'],second['s'])
            before=first['s']<start
            for k in ("rho","p","j","rad","drho","dp","dj","drad"):
                np.testing.assert_array_equal(first[k][before],both[k][before])
                np.testing.assert_array_equal(second[k][before],0.)
            increment=physical_difference(both,first)
            cross=physical_difference(increment,second)
            group={}
            for key,table in (("first",first),("second",second),("both",both),
                              ("added",increment),("interference",cross)):
                group[key]=dict(late=sector(table,start,READOUT),
                                full=sector(table,first["s"][0],READOUT))
            data[scheme]=group
        result[name]=data
    return dict(primary_window=ORIGINAL,late_window=LATE,readout=READOUT,
                units="Proper energy per normalized comoving volume, coefficient of y^2",
                cases=result)

def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo",type=Path,required=True)
    cli.add_argument("--tables",type=Path,required=True)
    cli.add_argument("--output",type=Path)
    cli.add_argument("--jobs",type=int,default=2)
    cli.add_argument("--steps",type=int,default=4096)
    cli.add_argument("--cutoff",type=int,choices=(120,240,480),default=120)
    cli.add_argument("--cases",nargs="+",choices=[g+"_"+a for g,a in CASES])
    cli.add_argument("--models",nargs="+",choices=MODELS)
    cli.add_argument("--read-only",action="store_true")
    args=cli.parse_args()
    cases=[tuple(c.split("_")) for c in args.cases] if args.cases else CASES
    if not args.read_only:
        generate(args.repo,args.tables,args.jobs,args.steps,cases,args.models,args.cutoff)
    if args.output:
        obj=summarize(args.tables,cases)
        obj["cutoff"]=args.cutoff
        args.output.write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")
if __name__=="__main__":main()
