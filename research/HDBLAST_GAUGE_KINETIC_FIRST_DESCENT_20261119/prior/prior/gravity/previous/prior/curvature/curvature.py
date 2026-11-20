"""Competing radiation equality and turnaround in the four-dimensional decay envelope."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
BRANCHES = ("0.5", "-0.5")


def evolve(rate, b0=0., r0=0., *, flat=False, tol=3e-12, method="DOP853"):
    """Positive curvature has b=Ka/M, rate=Gamma M/K^(3/2). Flat uses x=a(Gamma²/M)^(1/3)."""
    k = 0. if flat else 1.
    if min(rate, b0, r0)<0 or b0*k>=1+r0 or r0 >= 1:
        raise ValueError("An initially matter-dominated, expanding envelope is required")
    initial = [b0,math.sqrt(b0*(1+r0-k*b0)),1.,r0*b0,0.]

    def rhs(eta,y):
        b,v,p,q,u=y
        return (v, .5*p-k*b, -rate*b*p, rate*b*b*p, rate*b)

    def equality(eta,y):
        b,v,p,q,u=y
        return q/b-p if b else -1.

    def maximum(eta,y):
        return y[1]

    maximum.direction = -1
    maximum.terminal = True
    equality.direction = 1
    equality.terminal = flat
    sol=solve_ivp(rhs,(0.,12. if flat else math.pi+.02),initial,events=(equality,maximum),
                  method=method,rtol=tol,atol=tol*.002,max_step=.10)
    if not sol.success:
        raise RuntimeError(sol.message)
    def record(event):
        b,v,p,q,u=event
        return dict(b=float(b),velocity=float(v),gamma_elapsed=float(u),radiation_over_matter=float(q/(b*p)),
                    proper_acceleration=float(-p/(2*b*b)-q/b**3), H_over_Gamma=float(v/(rate*b*b)),
                    friedmann_residual=float(v*v+k*b*b-p*b-q))
    eq=record(sol.y_events[0][0]) if len(sol.y_events[0]) else None
    turn=record(sol.y_events[1][0]) if len(sol.y_events[1]) else None
    if not flat and turn is None: raise RuntimeError("First maximum not reached")
    if flat and eq is None: raise RuntimeError("Flat equality not reached")
    b,v,p,q,u=sol.y
    return dict(equality=eq,turn=turn,max_friedmann_residual=float(max(abs(v*v+k*b*b-p*b-q))))


def at_turn_by_angle(rate,b0=0.,r0=0.,tol=3e-12):
    """Independent evolution using the Friedmann constraint, arccot(aH/sqrt(K)) and the radiation fraction."""
    f0 = r0/(1+r0)
    start = math.asin(math.sqrt(b0/(1+r0)))
    def rhs(angle,state):
        f,u=state
        s=math.sin(angle)
        cot=math.cos(angle)/s if s else 0.
        b=math.exp(-u)*s*s/(1-f)
        return (2*(1-f)*(rate*b-cot*f)/(1+f),2*rate*b/(1+f))
    sol=solve_ivp(rhs,(start,math.pi/2),(f0,0.),method="DOP853",rtol=tol,atol=tol*.002,max_step=.06)
    if not sol.success: raise RuntimeError(sol.message)
    f,u=sol.y[:,-1]
    return float(f/(1-f)),float(u)


def critical(b0=0.,r0=0.,tol=3e-12,method="DOP853"):
    def residual(rate):
        return evolve(rate,b0,r0,tol=tol,method=method)["turn"]["radiation_over_matter"]-1
    return brentq(residual,.03,4.,xtol=3e-12,rtol=3e-12)


def universal():
    limit=critical()
    alternative=brentq(lambda g:at_turn_by_angle(g)[0]-1,.03,4.,xtol=3e-12)
    repeats={name:critical(tol=tol,method=method) for name,tol,method in (
        ("looser_dop853",1e-9,"DOP853"),("tight_rk45",2e-11,"RK45"))}
    samples=[]
    for g in sorted(set([.1,.25,.5,.75,1.,1.25,1.5,2.,3.,5.,8.,12.,16.,limit*(1-1e-5),limit*(1+1e-5)])):
        d=evolve(g)
        alt,u=at_turn_by_angle(g)
        samples.append(dict(g=g,ratio_at_turn=d["turn"]["radiation_over_matter"],
                        angle_check=alt,turn_gamma_elapsed=d["turn"]["gamma_elapsed"],
                        angle_gamma_elapsed=u,first_event="equality" if d["equality"] else "turnaround"))
    u_hi=brentq(lambda u:math.exp(u)-1-2*u,.1,2)
    flat=evolve(1.,flat=True)
    return dict(g_critical=limit,g_from_angle=alternative,refinements=repeats,
        k_over_GammaM_to_2thirds=limit**(-2/3),threshold_state=evolve(limit)["turn"],
        matter_start_sensitivity=[dict(b0=b,g_critical=critical(b)) for b in (1e-2,1e-4,1e-6)],
        scanned=samples,analytic_bounds=dict(g_no_equality_below=2*math.log(2)/math.pi,
            g_equality_guaranteed_above=2*u_hi*math.exp(u_hi)/(math.pi/2-1),
            gamma_elapsed_at_equality_lower=math.log(2),gamma_elapsed_at_equality_upper=u_hi),
        flat_early_matter=dict(**flat["equality"],h_over_Gamma=flat["equality"]["velocity"]/flat["equality"]["b"]**2),
        maximum_constraint_residual=max(evolve(g)["max_friedmann_residual"] for g in (.1,limit,3.)))


def anchor_scan(a,scalar,radiation,gamma):
    """Fix the material initial data and replace the spatial curvature, recomputing H_i by Friedmann."""
    r0=radiation/scalar
    hm=math.sqrt(scalar)
    M=a**3*scalar
    def b0(g):
        return (gamma/(hm*g))**(2/3)
    def residual(g, independent=False):
        if independent:
            return at_turn_by_angle(g,b0(g),r0)[0]-1.
        return evolve(g,b0(g),r0)["turn"]["radiation_over_matter"]-1.
    g=brentq(residual,.1,2.,xtol=4e-12)
    g_alt=brentq(lambda x:residual(x,True),.1,2.,xtol=4e-12)
    K=b0(g)*a*a*scalar
    flat0=b0(1.)
    feq=evolve(1.,flat0,r0,flat=True)["equality"]
    trials={}
    for mult in (.9,1.1):
        gg=g/mult**1.5
        trial=evolve(gg,b0(gg),r0)
        trials[str(mult)]=dict(first_event="radiation equality" if trial["equality"] else "turnaround",
            ratio_at_turn=trial["turn"]["radiation_over_matter"],
            H_over_Gamma_at_equality=trial["equality"]["H_over_Gamma"] if trial["equality"] else None)
    return dict(initial_scale=a,scalar_density=scalar,radiation_density=radiation,H_m=hm,
        H_i_flat=math.sqrt(scalar+radiation),Gamma_over_Htop=gamma,
        critical_g=g,critical_g_angle=g_alt,critical_K=K,
        maximum_allowed_K=a*a*(scalar+radiation),
        critical_fraction_of_initial_density=K/(a*a*(scalar+radiation)),curvature_controls=trials,
        flat_equality=dict(gamma_elapsed=feq["gamma_elapsed"],H_over_Gamma=feq["H_over_Gamma"],
                           a_over_start=feq["b"]/flat0,elapsed_in_Htop_inverse=feq["gamma_elapsed"]/gamma))


def analyze(background,widths):
    results={}
    for key in BRANCHES:
        name=f"phi_m={key}"
        old=background["branches"][name]
        width=widths["branches"][name]
        scales = {}
        for detuning in ("0.001","1e-10","1e-20"):
            point=width["scale_points"]["t="+detuning]
            gamma=point["Gamma_minimal_real_scalar_over_Htop"]
            microscopic=width["alpha"]**2*width["m_over_Htop"]**3*point["Htop_over_Mpl"]**2/(32*math.pi)
            if not math.isclose(gamma,microscopic,rel_tol=3e-14):
                raise ValueError("Inconsistent microscopic width normalization")
            scales[detuning]=anchor_scan(1.,(width["m_over_Htop"]/50)**2,0.,gamma)
        anchor=old["flat"]["harmonic"]["flat_anchor"]
        closed={}
        for model in ("off","always","harmonic"):
            trajectory=old["closed"][model]
            turn=trajectory["turnaround"]
            if turn is None or turn["radiation_fraction"]>=.5 or trajectory["first_harmonic_cycle"] is not None:
                raise ValueError("Closed-launch provenance inconsistent with measured history")
            closed[model]=dict(first_maximum_x=turn["x"],a_max=turn["a"],H_at_max=turn["h"],
                radiation_fraction=turn["radiation_fraction"],
                entire_expansion_radiation_fraction_upper=turn["radiation_fraction"]*(turn["a"]/trajectory["seed"]["a"])**4,
                friedmann_residual_at_turn=turn["scalar"]+turn["radiation"]-1./turn["a"]**2,
                Hprime_at_turn=1./turn["a"]**2-.5*turn["v"]**2-2*turn["radiation"],
                launch_K=1., launch_a=trajectory["seed"]["a"],
                launch_curvature_fraction=1./trajectory["seed"]["a"]**2/(trajectory["seed"]["scalar"]+trajectory["seed"]["radiation"]),
                completed_cycles=trajectory["completed_cycles"])
        results[name]=dict(m_over_Htop=width["m_over_Htop"],alpha=width["alpha"],
            flat_m_over_50=scales,measured_flat_anchor_at_t_0p001=anchor_scan(anchor["a"],anchor["scalar"],
            anchor["radiation"],old["gamma"]),closed_archived_launch_at_t_0p001=closed)
    return dict(model="4D dust decaying at constant rest width into free radiation",universal=universal(), branches=results)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background",type=Path,default=HERE.parent/"background/BACKGROUND.json")
    parser.add_argument("--widths",type=Path,default=HERE.parent/"widths/WIDTHS.json")
    parser.add_argument("--output",type=Path,default=HERE/"CURVATURE.json")
    args=parser.parse_args()
    result=analyze(json.loads(args.background.read_bytes()),json.loads(args.widths.read_bytes()))
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")


if __name__=="__main__":
    main()
