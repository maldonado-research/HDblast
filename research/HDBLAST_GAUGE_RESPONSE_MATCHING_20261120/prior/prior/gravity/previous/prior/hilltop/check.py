"""Check the closed hilltop transfer with its exact linear gain and full nonlinear launch."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import hyp2f1


def engine(path):
    spec=importlib.util.spec_from_file_location("hdblast_hilltop_reference",path/"hilltop.py")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def local_force(p,top):
    # Use the cubic on the outgoing side. The polynomial for U' has an exactly
    # zero constant at the stationary point. Drop it before division by the offset.
    i=int(np.searchsorted(p.theta,top+p.sign*1e-6)-1)
    c=p.spline.c[:,i]
    z=Polynomial([top-p.theta[i],1.])
    phi,f=(Polynomial(cc[::-1])(z) for cc in c.T)
    P=(1-phi/p.pm)**2*(1+p.b*phi)
    Q=f*f*p.v0
    D=P.deriv()*Q-P*Q.deriv()
    reduced=Polynomial(D.coef[1:])
    def force(offset):
        return float(reduced(offset)/Q(offset)**2)
    return force,float(D(0.)/Q(0.)**2)


def direct_to(p,hill,logseed,offset_end,rtol):
    top=hill.top
    force,root_residual=local_force(p,top)
    transition=2e-5
    at_switch=p.sign*transition
    mismatch=force(at_switch)-p.at(top+at_switch)[1]/at_switch
    if abs(mismatch)>1e-6:raise ValueError("Hilltop derivative charts disagree")
    def rhs(x,y):
        l,g,loga,h=y
        d=p.sign*math.exp(l)
        f=force(d) if abs(d)<transition else p.at(top+d)[1]/d
        return (g,-3*h*g-g*g-3*f,h,math.exp(-2*loga)-.5*(d*g)**2)
    def stop(x,y):return y[0]-math.log(offset_end)
    stop.terminal=True;stop.direction=1
    d0=p.sign*math.exp(logseed)
    state=[logseed,0.,-.5*math.log(p.at(top+d0)[0]),0.]
    sol=solve_ivp(rhs,(0.,50.),state,method="DOP853",max_step=.04,rtol=rtol,atol=rtol*.002,
                  events=stop)
    if not sol.success or len(sol.t_events[0])!=1:raise ValueError("Nonlinear departure missing")
    x=float(sol.t_events[0][0])
    l,g,loga,h=sol.y_events[0][0]
    delta=p.sign*offset_end
    s=np.array([top+delta,delta*g,loga,0.])
    Hc=math.sqrt(p.at(s[0])[0]+s[1]*s[1]/6-math.exp(-2*loga))
    return x,s,dict(friedmann_at_departure=float(h-Hc), force_chart_difference=mismatch,
                     stationary_residual_in_derivative=root_residual)


def matched_to(p,hill,logseed,offset_end):
    x=brentq(lambda t:hill.growth(t)[0]-(math.log(1e-5)-logseed),0.,41.)
    start,_=hill.linear_state(x,1e-5)
    def rhs(t,y):
        th,v,loga=y
        u,du,_=p.at(th)
        h=math.sqrt(u+v*v/6-math.exp(-2*loga))
        return [v,-3*h*v-3*du,h]
    def stop(t,y):return p.sign*(y[0]-hill.top)-offset_end
    stop.terminal=True;stop.direction=1.
    sol=solve_ivp(rhs,(x,x+2.),start[:3],rtol=2e-12,atol=1e-15,method="DOP853",
                  events=stop,max_step=.017)
    if not sol.success or len(sol.t_events[0])!=1:raise ValueError("Matched departure missing")
    return float(sol.t_events[0][0]),sol.y_events[0][0]


def linear_closed(hill,x):
    # With z=−sinh²(hx), the regular even-in-time solution is Gauss hypergeometric.
    disc=math.sqrt(9+4*hill.mu/hill.h**2)
    a,b=(3+disc)/4,(3-disc)/4
    z=-math.sinh(hill.h*x)**2
    F=hyp2f1(a,b,.5,z)
    dFdz=2*a*b*hyp2f1(a+1,b+1,1.5,z)
    g=dFdz/F*(-hill.h*math.sinh(2*hill.h*x))
    return np.array([math.log(F),g])


def calculate(study,archive,history):
    M=engine(study)
    mod,ref=M.BG.source(archive)
    tables=mod.build_tables()
    result={}
    for key,branch in history["branches"].items():
        p=M.BG.Potential(mod,tables,float(key.split("=")[1]),ref[key])
        hill=M.ClosedLaunch(p)
        dets={}
        for det,run in branch["detunings"].items():
            primary=run["primary"]
            x=primary["x_match"]
            logseed=primary["log_abs_bounce_seed"]
            linear_error=linear_closed(hill,x)-hill.growth(x)
            # 0.03 Mpl is far beyond the 10^-5 Mpl linear-to-nonlinear overlap.
            t,s,diagnostics=direct_to(p,hill,logseed,.03,2e-12)
            tm,sm=matched_to(p,hill,logseed,.03)
            row=dict(linear_log_gain_error=float(linear_error[0]),
                     linear_log_derivative_error=float(linear_error[1]),
                     departure_abs_field=.03, nonlinear_departure_x=t,
                     matched_minus_nonlinear_x=tm-t,
                     matched_minus_nonlinear_loga=float(sm[2]-s[2]),
                     matched_relative_velocity=float(sm[1]/s[1]-1),**diagnostics)
            if max(abs(linear_error))>1e-8 or max(abs(tm-t),abs(sm[2]-s[2]),abs(sm[1]/s[1]-1))>1e-5:
                raise ValueError("Independent nonlinear departure disagrees")
            row["absolute_seed_in_Mpl"]=float(math.exp(logseed))
            row["slope_shift_matching_seed"]=float(hill.mu*math.exp(logseed)/3)
            if float(det)==.001:
                def margin(logseed):
                    tx,initial,_=direct_to(p,hill,logseed,.03,2e-12)
                    case=M.ClosedLaunch(p)
                    case.linear_state=lambda x0,match:(initial,logseed)
                    return case.anchors(tx,.03,run["gamma"],rtol=2e-11,short=True)["anchors"]["24_cross"]["log_margin"]
                row["actual_bounce_log_margin_for_primary_seed"]=margin(logseed)
                b=run["brackets_all_handoff_points"]
                row["actual_bounce_bracket_margins"]=[margin(float(item["log10_seed_over_archive"])*math.log(10.)+
                                                             math.log(run["archive_displacement_scale"])) for item in b]
                m0,m1=row["actual_bounce_bracket_margins"]
                if m0>=0 or m1<=0:raise ValueError("Nonlinear physical bounce lost the threshold bracket")
            dets[det]=row
        result[key]=dets
    return result

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive",type=Path,required=True)
    parser.add_argument("--study",type=Path,default=Path(__file__).resolve().parent,
                        help="Directory containing hilltop.py and HILLTOP.json")
    parser.add_argument("--input",type=Path,help="Override the registered hilltop JSON")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    history=json.loads((args.input or args.study/"HILLTOP.json").read_text())
    args.output.write_text(json.dumps(calculate(args.study,args.archive,history),indent=2,allow_nan=False)+"\n")
