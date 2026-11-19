"""Can the closed hilltop bounce supply enough expansion for the calibrated decay?"""
import argparse
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


HERE=Path(__file__).resolve().parent.parent


def load(name):
    path=HERE/name/(name+".py")
    spec=importlib.util.spec_from_file_location("modulus_"+name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BG=load("background")
ENVELOPE=load("curvature")


def field_top(pot):
    top=brentq(lambda theta:pot.at(theta)[1],-1e-3,1e-3,xtol=8e-16)
    # Differentiate the spline analytically at its stationary point.
    phi,f=pot.spline(top)
    p1,f1=pot.spline(top,1)
    p2,f2=pot.spline(top,2)
    z=1-phi/pot.pm
    j=1+pot.b*phi
    first=-2/(pot.pm*z)+pot.b/j
    second=-2/(pot.pm*z)**2-(pot.b/j)**2
    u=pot.at(top)[0]
    logcurv=second*p1*p1+first*p2-2*(f2/f-(f1/f)**2)
    curvature=u*logcurv
    h=math.sqrt(u)
    mu=-3*curvature
    def riccati(x,state):
        growth=state[1]
        return growth,mu-3*h*math.tanh(h*x)*growth-growth*growth
    sol=solve_ivp(riccati,(0,42.),[0.,0.],rtol=2e-12,atol=2e-14,
                  method="DOP853",max_step=.06,dense_output=True)
    if not sol.success or mu<=0:
        raise ValueError("Closed hilltop linearization")
    return top,h,mu,sol.sol


def constrained(pot,state,t0,t1,gamma,rtol,turning=True):
    """On the expanding portion, solve Friedmann's constraint at every RHS evaluation."""
    def equation(t,y):
        th,v,loga,r=y
        u,du,_=pot.at(th)
        h=math.sqrt(max(0.,u+v*v/6+r-math.exp(-2*loga)))
        return [v,-(3*h+gamma)*v-3*du,h,gamma*v*v/3-4*h*r]
    def crossing(t,y):return pot.sign*(y[0]-pot.tmin)
    def turn(t,y):return y[1]
    crossing.direction=1.
    events=(crossing,turn) if turning else None
    result=solve_ivp(equation,(t0,t1),state,method="DOP853",rtol=rtol,
                     atol=rtol*np.array([.01,.01,.1,1e-6]),events=events,
                     max_step=2*math.pi/pot.mharm/12,dense_output=True)
    if not result.success:raise ValueError(result.message)
    th,v,loga,r=result.y
    u=np.array([pot.at(value)[0] for value in th])
    if np.min(u+v**2/6+r-np.exp(-2*loga))<=0:
        raise ValueError("Nonlinear evolution touched a non-expanding branch before matching")
    return result


def values(pot,y):
    th,v,loga,r=y
    u=pot.at(th)[0]
    scalar=u+v*v/6
    total=scalar+r-math.exp(-2*loga)
    if total<=0:raise ValueError("Handoff beyond first turnaround")
    return dict(theta=float(th),v=float(v),loga=float(loga),a=math.exp(loga),
                scalar=float(scalar),radiation=float(r),h=math.sqrt(total))


def gate_time(pot,segment):
    history=[]
    for start,stop in zip(segment.t_events[0],segment.t_events[0][1:]):
        turns=[float(v) for v in segment.t_events[1] if start<v<stop]
        if len(turns)!=2:raise ValueError("Cycle lost a turning point")
        t=np.sort(np.r_[np.linspace(start,stop,129),turns])
        states=segment.sol(t)
        h=[]; errors=[]
        for y in states.T:
            res=values(pot,y); h.append(res["h"]/pot.mharm)
            dx=y[0]-pot.tmin
            if abs(dx)>1e-8:
                errors.append(abs(6*pot.at(y[0])[0]/(pot.mharm*dx)**2-1))
        history.append(max(h)<=.01 and max(errors)<=.1)
        if len(history)>1 and history[-1] and history[-2]:return float(stop)
    raise ValueError("No two qualifying harmonic cycles")


def critical_rate(b,r):
    def residual(rate):return ENVELOPE.at_turn_by_angle(rate,b,r,tol=5e-11)[0]-1
    return brentq(residual,.2,2.,xtol=4e-10)


class ClosedLaunch:
    def __init__(self,pot):
        self.p=pot
        self.top,self.h,self.mu,self.growth=field_top(pot)

    def linear_state(self,x,match):
        p=self.p
        loggain,z=self.growth(x)
        offset=p.sign*match
        logseed=math.log(match)-loggain
        loga=math.log(math.cosh(self.h*x)/self.h)
        return np.array([self.top+offset,offset*z,loga,0.]),float(logseed)

    def anchors(self,x,match,gamma,rtol=2e-10,short=False):
        p=self.p
        start,logseed=self.linear_state(x,match)
        init=constrained(p,start,x,x+18,0.,rtol)
        first=float(init.t_events[0][0])
        stage=gate_time(p,init)
        data=constrained(p,init.sol(stage),stage,first+(27 if short else 40),gamma,rtol)
        points={}
        for age in ([24] if short else [12,24,36]):
            cross=next(t for t in data.t_events[0] if t>=first+age)
            if short:
                samples=[("cross",cross)]
            else:
                turn=next(t for t in data.t_events[1] if t>cross)
                samples=[("cross",cross),("turn",turn)]
            for name,t in samples:
                pt=values(p,data.sol(t));a=pt["a"];e=pt["scalar"];r=pt["radiation"]/e
                g=gamma*a**3*e; b=1/(a*a*e); crit=critical_rate(b,r)
                if not b<1+r:raise ValueError("Invalid envelope start")
                points[f"{age}_{name}"]=dict(**pt,x=float(t),g=g,b=b,r=r,
                                               critical_rate=crit,log_margin=math.log(g/crit))
        return dict(x_match=x,log_abs_bounce_seed=logseed,first_min_x=first,width_switch_x=stage,
                    loga_at_first_min=values(p,init.sol(first))["loga"],anchors=points)

    def overlap(self):
        p=self.p; offset=p.sign*2e-7
        m=2e-4
        x=brentq(lambda x:self.growth(x)[0]-math.log(m/abs(offset)),0.,8.)
        matched,logseed=self.linear_state(x,m)
        # Full time-symmetric bounce from constraint, with h evolved through zero.
        th0=self.top+offset;u0=p.at(th0)[0]
        init=np.array([th0,0.,-.5*math.log(u0),0.,0.,0.,0.])
        sol=BG.segment(p,1,0.,init,0.,x+.05,1e-12)
        hit=brentq(lambda t:p.sign*(sol.sol(t)[0]-self.top)-m,x-.015,x+.015)
        native=sol.sol(hit)
        return dict(checked_abs_seed=abs(offset),comparison_abs_field=m,linear_event=x,
                    direct_event=float(hit),match_time_error=x-hit,relative_velocity_error=float(matched[1]/native[1]-1),
                    loga_error=float(matched[2]-native[2]),reconstructed_seed=math.exp(logseed))


def calculate_branch(p,dets,controls):
    case=ClosedLaunch(p)
    records={}
    for det in dets:
        flags=controls and det==.001
        gamma=p.gamma*det/.001
        sigma=p.dq*math.sqrt(det/.001)
        guess=math.log(1/gamma)/3+0.2
        def threshold(match):
            @lru_cache(maxsize=None)
            def trial(x):return case.anchors(x,match,gamma,short=True)["anchors"]["24_cross"]["log_margin"]
            lo,hi=guess-1.6,guess+2.2
            # Move towards equality from well-defined postinflationary expanding histories.
            while trial(lo)>0:lo-=.35
            while trial(hi)<0:hi+=.35
            x=brentq(trial,lo,hi,xtol=1.5e-6)
            y=case.anchors(x,match,gamma)
            slope=(trial(x+.015)-trial(x-.015))/.030
            log_logseed_slope=-case.growth(x)[1]/math.log(10.)
            h=y["anchors"]["24_cross"]
            # Local root estimate for alternative matching ages. Hard brackets checked below.
            shifts={name:-row["log_margin"]/slope*log_logseed_slope for name,row in y["anchors"].items()}
            return dict(**y,log10_bounce_seed_over_archive=(y["log_abs_bounce_seed"]-math.log(sigma))/math.log(10.),
                        dlog_margin_dx =slope,estimated_log10_seed_shifts=shifts,
                        angle_ratio_at_turn=ENVELOPE.at_turn_by_angle(h["g"],h["b"],h["r"])[0],
                        direct_ratio_at_turn=ENVELOPE.evolve(h["g"],h["b"],h["r"])["turn"]["radiation_over_matter"])
        res=threshold(1e-5)
        root=res["x_match"]
        bracket=[]
        for x in [root-.0018,root+.0018]:
            hit=case.anchors(x,1e-5,gamma)
            margins=[v["log_margin"] for v in hit["anchors"].values()]
            h=hit["anchors"]["24_cross"]
            direct=ENVELOPE.evolve(h["g"],h["b"],h["r"])
            bracket.append(dict(first_event="radiation equality" if direct["equality"] else "turnaround",
                                direct_ratio_at_first_turn=direct["turn"]["radiation_over_matter"],
                                log10_seed_over_archive=(hit["log_abs_bounce_seed"]-math.log(sigma))/math.log(10.),
                                min_log_margin=min(margins),max_log_margin=max(margins)))
        if (bracket[0]["max_log_margin"]>=0 or bracket[1]["min_log_margin"]<=0
                or bracket[0]["first_event"]!="turnaround" or bracket[1]["first_event"]!="radiation equality"):
            raise ValueError("Handoff ages don't share a threshold bracket")
        other={}
        if flags:
            for m in [1e-6,1e-4]:
                alt=threshold(m)
                other[str(m)]=dict(x_match=alt["x_match"],log10_seed_over_archive=alt["log10_bounce_seed_over_archive"])
        records[f"{det:g}"]=dict(gamma=gamma,archive_displacement_scale=sigma,primary=res,
                                 brackets_all_handoff_points=bracket,departure_control=other)
    growth=(-3*case.h+math.sqrt(9*case.h**2+4*case.mu))/2
    predicted=growth/(3*case.h)-.5
    scaling=dict(predicted_pow_t_for_relative_seed=predicted)
    if "1e-10" in records and "1e-20" in records:
        a,b=(records[key]["primary"]["log10_bounce_seed_over_archive"] for key in ("1e-10","1e-20"))
        measured=(a-b)/10
        if abs(measured-predicted)>1e-4:
            raise ValueError("Low-width asymptotic failure")
        scaling["measured_pow_t_for_relative_seed"]=measured
    return dict(stationary_theta=case.top,negative_mass_squared_over_Htop_squared=case.mu,
                limiting_linear_growth_rate=growth,scaling=scaling,overlap=case.overlap(),detunings=records)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--detunings",nargs="+",type=float,default=[.001,1e-10,1e-20])
    parser.add_argument("--minimum",nargs="+",type=float,default=[.5,-.5])
    parser.add_argument("--quick",action="store_true")
    args=parser.parse_args()
    mod,ref=BG.source(args.archive)
    tables=mod.build_tables()
    result={}
    for pm in args.minimum:
        p=BG.Potential(mod,tables,pm,ref[f"phi_m={pm:g}"])
        result[f"phi_m={pm:g}"]=calculate_branch(p,args.detunings,not args.quick)
        print(f"Completed hilltop {pm:g}",flush=True)
    args.output.write_text(json.dumps(dict(branches=result),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":main()
