"""Two-endpoint variance on four early compensated closed bounces at t=10^-3."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def propagate(p, top, engine, velocity, center, rtol, step):
    wave = 8  # n=2
    theta0 = top+center
    u, du, _ = p.at(theta0)
    a0 = 1/math.sqrt(u+velocity**2/6)
    alpha0 = 2*(wave-3)/a0**2
    chart = 1-velocity**2/alpha0
    if chart <= 0:
        raise ValueError("The specified (q,q') chart must have positive orientation")
    D0 = alpha0*chart
    kappa = p.Hpl**2/(2*math.pi**2*a0**3)
    kappa_eff = kappa*chart
    # Two Newtonian columns in (q0,q0') coordinates and two others in (q0,P0).
    # The canonical pair's P=q'-vPhi needs only alpha>0 for reconstruction.
    q0, w0 = np.eye(2)
    phi0 = (-3*du*q0-velocity*w0)/D0
    init = np.r_[[theta0, velocity, math.log(a0), 0.], q0,w0,phi0,q0,w0]
    def quantities(y):
        theta,v,l,h = y[:4]
        u,du,_ = p.at(theta)
        inv2 = math.exp(-2*l)
        vp = -3*(h*v+du)
        alpha = 2*(wave-3)*inv2
        omega = wave*inv2+engine.local_mass(p,theta)
        return u,du,inv2,vp,alpha,omega
    def equation(t,y):
        _,v,_,h = y[:4]
        _,du,inv2,vp,alpha,omega = quantities(y)
        q,w,phi,qc,P = y[4:].reshape(5,2)
        phip = -h*phi + .5*v*q
        wp = -3*h*w-omega*q+4*v*phip-6*du*phi
        phic = (vp*qc-v*P)/alpha
        return np.r_[[v,vp,h,inv2-.5*v*v],
                     w,wp,phip,P+v*phic,-3*h*P-(omega-1.5*v*v)*qc+vp*phic]
    def horizon(t,y):
        # At positive B, the unit-RMS nine-harmonic floor is sqrt(9 kappa_eff B/2).
        return 9*kappa_eff*y[5]/2-.01**2
    horizon.terminal = True
    horizon.direction = 1
    sol = solve_ivp(equation,(0.,1.),init,method="DOP853",rtol=rtol,atol=rtol*1e-3,
                    max_step=step,events=horizon,dense_output=True)
    if not sol.success or len(sol.t_events[0])!=1:
        raise RuntimeError("Two-endpoint floor did not reach .01")
    t = float(sol.t_events[0][0])
    def wedge(a,b):
        return a[0]*b[1]-a[1]*b[0]
    diagnostics = []
    for state in sol.sol(np.linspace(0,t,129)).T:
        _,v,l,h = state[:4]
        u,_,inv2,vp,alpha,_ = quantities(state)
        q,w,phi,qc,P = state[4:].reshape(5,2)
        ham_terms = np.stack(((alpha-v*v)*phi,v*w,-vp*q))
        ham = np.linalg.norm(ham_terms.sum(axis=0))/np.linalg.norm(abs(ham_terms).sum(axis=0))
        # Both initial charts map to the same physical observables.
        canonical_q = qc[0]*q0 + qc[1]*(w0-velocity*phi0)
        area_ratio = math.exp(3*l)/a0**3*chart*wedge(q,w-v*phi)
        diagnostics.append([ham,abs(h*h-u-v*v/6+inv2),abs(area_ratio-1),
                            np.linalg.norm(canonical_q-q)/np.linalg.norm(q),w[1]])
    checks = np.array(diagnostics)
    y = sol.y_events[0][0]
    b = float(y[5])
    canonical_b = float(y[11])
    if max(np.max(checks[:,:4],axis=0))>1e-7 or min(checks[:,4])<=0:
        raise ArithmeticError("Compensated velocity transfer audit failed")
    return dict(velocity=velocity,center_offset=center,bounce_radius=a0,kappa_in_q_P=kappa,
                q_qprime_chart_factor=chart,kappa_in_q_qprime=kappa_eff,horizon_x=t,
                homogeneous_offset_at_horizon=float(y[0]-top),B_in_q_qprime=b,B_in_q_P=canonical_b,
                rolling_rms_floor=math.sqrt(9*kappa_eff*b/2),
                canonical_rms_floor=math.sqrt(9*kappa*canonical_b/2),
                max_scaled_Hamiltonian_constraint=max(checks[:,0]),
                max_Friedmann_error=max(checks[:,1]),max_relative_symplectic_drift=max(checks[:,2]),
                max_relative_canonical_transfer_difference=max(checks[:,3]),
                minimum_sampled_B_prime=min(checks[:,4]))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prior",type=Path,required=True)
    cli.add_argument("--manifold",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    args = cli.parse_args()
    H = module(args.prior/"hilltop"/"hilltop.py","velocity_parent")
    rolling = module(Path(__file__).with_name("rolling.py"),"rolling_helper")
    archive, records = H.BG.source(args.prior/"archive")
    tables = archive.build_tables()
    inputs = json.loads(args.manifold.read_text())
    out = {}
    for key,branch in inputs["branches"].items():
        p = H.BG.Potential(archive,tables,float(key.split("=")[1]),records[key])
        top = H.field_top(p)[0]
        values = []
        for row in branch["samples"]:
            v, center = row["v"],row["center"]
            if not math.isclose(abs(v),math.sqrt(.06)):
                continue
            result = propagate(p,top,rolling,v,center,1e-11,.007)
            fine = propagate(p,top,rolling,v,center,2e-12,.003)
            # A resolvable perturbation larger than the reference width and center error.
            shifted = [propagate(p,top,rolling,v,center+delta,2e-12,.003)["horizon_x"]
                       for delta in [-1e-9,1e-9]]
            result.update(abs_center_bracket_shift=1e-9,
                          calibrated_distance_to_center=row["conditional_width"],
                          prior_paired_match_x=row["match_time"],
                          max_horizon_change_under_center_bracket=max(abs(x-fine["horizon_x"]) for x in shifted),
                          horizon_refinement_difference=abs(fine["horizon_x"]-result["horizon_x"]))
            values.append(result)
        if len(values)!=2:
            raise ValueError("Missing either extreme launch")
        out[key] = values
    args.output.write_text(json.dumps(dict(detuning=.001,rms_goal=.01,branches=out),
                                      indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
