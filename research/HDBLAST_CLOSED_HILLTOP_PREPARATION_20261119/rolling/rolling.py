"""Classical Newtonian-gauge transport for scalar harmonics of a tuned closed bounce."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def local_mass(p, theta):
    index = int(np.searchsorted(p.theta, theta, side="right") - 1)
    z = theta - p.theta[index]
    c3, c2, c1, c0 = p.coeff[:, index]
    phi, f = ((c3*z + c2)*z + c1)*z + c0
    dphi, df = (3*c3*z + 2*c2)*z + c1
    ddphi, ddf = 6*c3*z + 2*c2
    u = (1-phi/p.pm)**2 * (1+p.b*phi) / (f*f*p.v0)
    L = -2/(p.pm-phi) + p.b/(1+p.b*phi)
    dL = -2/(p.pm-phi)**2 - (p.b/(1+p.b*phi))**2
    slope = L*dphi - 2*df/f
    return 3*u*(slope*slope + dL*dphi*dphi + L*ddphi - 2*(ddf/f - (df/f)**2))


def integrate(p, hill, root_force, seed, ns, endpoints, rtol, step, detuning, match_x):
    """Run momentum and spatial-trace constrained solvers alongside the stationary field."""
    top, hs = hill.top, hill.h
    growth = .5*(math.sqrt(9*hs**2 + 4*hill.mu) - 3*hs)
    k2 = (ns*(ns+2))[:, None]
    count = len(ns)
    def force(d):
        return root_force(d) if abs(d) < 2e-5 else p.at(top+d)[1]/d
    d0 = p.sign * math.exp(seed)
    a0 = p.at(top+d0)[0]**-.5
    initial = np.zeros((11,count,2))
    # q, w=q'/B, P=Phi/(d B), B=exp(growth*x); q on the RHS is Q/B.
    initial[[0,3,7,9],:,0] = 1.
    initial[[1,4,8,10],:,1] = 1.
    initial[[2,5],:,0] = (-3*force(d0)*a0**2/(2*(k2[:,0]-3)))
    def bundle(x, state):
        l, g, loga, h = state[:4]
        d = p.sign * math.exp(l)
        f = force(d)
        hprime = math.exp(-2*loga) - .5*(g*d)**2
        return d, f, hprime, math.exp(-2*loga), local_mass(p, top+d)
    def equation(x, state):
        _, g, _, h = state[:4]
        d, f, hprime, ainv2, mass = bundle(x, state)
        q,w,P,qt,wt,Pt,R,qs,ws,qf,wf = state[4:].reshape(11,count,2)
        # Momentum solution sets Phi'/(dB) algebraically. Spatial trace evolves it.
        velocity_P = -h*P + .5*g*q
        friction = 3*h + growth
        omega = k2*ainv2 + mass
        def acceleration(q, w, P, R):
            return -friction*w - omega*q + d*d*(4*g*R - 6*f*P)
        stationary_a2inv = (hs/math.cosh(hs*x))**2
        rows = (w-growth*q, acceleration(q,w,P,velocity_P),
                velocity_P-(g+growth)*P,
                wt-growth*qt, acceleration(qt,wt,Pt,R), R-(g+growth)*Pt,
                -(4*h+g+growth)*R-(2*hprime+3*h*h-ainv2+.5*(d*g)**2)*Pt
                +.5*(g*wt-3*f*qt),
                ws-growth*qs,
                -(3*hs*math.tanh(hs*x)+growth)*ws - (k2*stationary_a2inv-hill.mu)*qs,
                wf-growth*qf, -friction*wf-omega*qf)
        return np.r_[[g,-3*h*g-g*g-3*f,h,hprime], np.asarray(rows).ravel()]
    def event(m, final):
        def hit(x, state):
            return state[0] - math.log(m)
        hit.terminal = final
        hit.direction = 1
        return hit
    sol = solve_ivp(equation, (0.,50.), np.r_[[seed,0.,math.log(a0),0.],initial.ravel()],
                    method="DOP853",rtol=rtol,atol=rtol*1e-4,max_step=step,
                    events=[event(m, m==endpoints[-1]) for m in endpoints], dense_output=True)
    if not sol.success or any(len(y)!=1 for y in sol.t_events):
        raise ValueError(sol.message)
    def relrow(values, ref):
        return np.linalg.norm(values-ref,axis=1)/np.maximum(np.linalg.norm(ref,axis=1),1e-250)
    def normalized_residual(x, state):
        _, g, loga, h = state[:4]
        d,f,_,ainv2,_ = bundle(x, state)
        q,w,P,qt,wt,Pt,R,_,_,_,_ = state[4:].reshape(11,count,2)
        A = (2*(k2-3)*ainv2 - (d*g)**2)*P
        B = g*w
        D = (3*h*g+3*f)*q
        scale = np.linalg.norm(abs(A)+abs(B)+abs(D),axis=1)
        ham = np.linalg.norm(A+B+D,axis=1)/np.maximum(scale,1e-250)
        trace_mom = R+h*Pt - .5*g*qt
        scale = np.linalg.norm(abs(R)+abs(h*Pt)+abs(.5*g*qt),axis=1)
        mom = np.linalg.norm(trace_mom,axis=1)/np.maximum(scale,1e-250)
        fried = abs(h*h - (d*g)**2/6 - p.at(top+d)[0] + ainv2)
        return float(max(ham)),float(max(mom)),float(fried)
    errors = np.max([normalized_residual(x,y) for x,y in zip(sol.t,sol.y.T)],axis=0)
    data = []
    for end, te, ye in zip(endpoints, sol.t_events,sol.y_events):
        x, state = float(te[0]), ye[0]
        _,g,loga,h = state[:4]
        d,f,_,ainv2,mass = bundle(x,state)
        q,w,P,qt,wt,Pt,_,qs,ws,qf,wf = state[4:].reshape(11,count,2)
        B = math.exp(growth*x)
        errors_q,errors_w,errors_p = relrow(qt,q),relrow(wt,w),relrow(Pt,P)
        dq, dw = relrow(q,qs),relrow(w,ws)
        metric_q, metric_w = relrow(q,qf), relrow(w,wf)
        out = []
        for i,n in enumerate(ns):
            out.append(dict(n=int(n),delta_theta=(B*q[i]).tolist(),delta_theta_prime=(B*w[i]).tolist(),
                            phi=(B*d*P[i]).tolist(),stationary_delta_theta=(B*qs[i]).tolist(),
                            stationary_delta_theta_prime=(B*ws[i]).tolist(),
                            fractional_change_delta_theta=float(dq[i]),fractional_change_delta_theta_prime=float(dw[i]),
                            metric_source_fractional_change_delta_theta=float(metric_q[i]),
                            metric_source_fractional_change_delta_theta_prime=float(metric_w[i]),
                            phi_row_norm_over_delta_theta_row_norm=float(abs(d)*np.linalg.norm(P[i])/np.linalg.norm(q[i])),
                            trace_vs_momentum_fractional_difference=float(max(errors_q[i], errors_w[i], errors_p[i]))))
        data.append(dict(abs_offset=end,x=x,loga=loga,h=h,v=d*g,mass_squared_over_Href_squared=mass,
                         max_delta_theta_row_change=max(dq),max_delta_theta_prime_row_change=max(dw),
                         modes=out))
    # Newtonian-gauge ADM symplectic form per real unit-RMS harmonic at the bounce.
    # Degeneracy nine enters the whole n=2 multiplet RMS, not its mode commutator.
    href = p.Hpl*math.sqrt(detuning/1e-3)
    kappa = href**2 / (2*math.pi**2*a0**3)
    stationary_kappa = href**2 * hs**3 / (2*math.pi**2)
    def snapshot(x):
        state = sol.sol(x)
        rows = state[4:].reshape(11,count,2)
        factor = math.exp(growth*x)
        b = factor*rows[0,0,1]
        b_stationary = factor*rows[7,0,1]
        return dict(x=x, abs_background_offset=math.exp(state[0]), B_rolling=b,
                    B_stationary=b_stationary,
                    n2_multiplet_rms_floor=math.sqrt(9*kappa*abs(b)/2),
                    stationary_n2_multiplet_rms_floor=math.sqrt(9*stationary_kappa*abs(b_stationary)/2),
                    single_real_harmonic_rms_floor=math.sqrt(kappa*abs(b)/2))
    final = sol.t[-1]
    rms_goal = .01
    def horizon(key):
        return brentq(lambda x:snapshot(x)[key]-rms_goal,0.,final,xtol=2e-14)
    rolling_horizon = horizon("n2_multiplet_rms_floor")
    stationary_horizon = horizon("stationary_n2_multiplet_rms_floor")
    uncertainty = dict(benchmark_rms=rms_goal, bounce_radius=a0,
                       Href_over_Mpl=href,kappa_per_real_harmonic=kappa,
                       rolling_horizon_x=rolling_horizon,stationary_horizon_x=stationary_horizon,
                       samples={name:snapshot(x) for name,x in
                                [("below",.99*stationary_horizon),("at_stationary_horizon",stationary_horizon),
                                 ("above",1.01*stationary_horizon),("at_rolling_horizon",rolling_horizon),
                                 ("at_prior_match",match_x)]})
    audits = []
    for n in (2,8):
        if n > ns[-1]:
            continue
        wave = n*(n+2)
        def ratio(x):
            l, g, loga, _ = sol.sol(x)[:4]
            return (math.exp(l+loga)*g)**2 / (2*(wave-3))
        hi = min(2.5,ratio(final))
        def ratio_time(value):
            return brentq(lambda x:ratio(x)-value,0.,final,xtol=2e-14)
        x0 = ratio_time(min(.4,hi/4))
        x1 = final if hi == ratio(final) else ratio_time(hi)
        start = sol.sol(x0)
        l, g, loga0, h = start[:4]
        d, f, _, ainv2, _ = bundle(x0, start)
        v, vp, alpha = d*g, -3*d*(h*g+f), 2*(wave-3)*ainv2
        # Scale Phi by 1/sqrt(alpha_start) and choose a balanced basis for the
        # two-dimensional plane of constrained initial data even at D=0.
        phi_scale = 1/math.sqrt(alpha)
        normal = np.array([-vp,v,(alpha-v*v)*phi_scale])
        init = np.linalg.svd(normal.reshape(1,3),full_matrices=True)[2][1:].T
        def local_equation(x, state):
            bg = sol.sol(x)
            _, g, _, h = bg[:4]
            d, f, _, ainv2, mass = bundle(x,bg)
            q, w, z = state.reshape(3,2)
            zdot = -h*z + .5*d*g*q/phi_scale
            wdot = -3*h*w-(wave*ainv2+mass)*q + d*phi_scale*(4*g*zdot-6*f*z)
            return np.asarray([w,wdot,zdot]).ravel()
        local = solve_ivp(local_equation,(x0,x1),init.ravel(),method="DOP853",
                          max_step=(x1-x0)/12,rtol=rtol,atol=rtol*1e-4,dense_output=True)
        if not local.success:
            raise ValueError(local.message)
        def wedge(a, b):
            return a[0]*b[1]-a[1]*b[0]
        def invariants(x):
            l,g,loga,_ = sol.sol(x)[:4]
            q,w,z = local.sol(x).reshape(3,2)
            volume = math.exp(3*(loga-loga0))
            plain = volume*wedge(q,w)
            corrected = plain - volume*(p.sign*math.exp(l)*g)*phi_scale*wedge(q,z)
            return corrected,plain
        initial_area,plain0 = invariants(x0)
        grid = np.array([invariants(t) for t in np.linspace(x0,x1,65)])
        err = max(abs(grid[:,0]/initial_area-1))
        audits.append(dict(n=n,x_start=x0,x_stop=x1,ratio_at_start=ratio(x0),ratio_at_stop=ratio(x1),
                           symplectic_max_relative_drift=float(err),
                           field_only_volume_wronskian_change_over_initial_symplectic_area=
                           float(max(abs(grid[:,1]-plain0))/abs(initial_area))))
        if err>1e-7:
            raise ArithmeticError("Constraint-surface symplectic check failed")
    return dict(minus_stationary_mass_squared_over_Href_squared=hill.mu,stationary_hubble_scale_squared_over_Href_squared=hs*hs,
                abs_seed=abs(d0),constraint_max_scaled_hamiltonian=errors[0],
                independent_trace_max_scaled_momentum=errors[1],max_friedmann_error=errors[2],endpoints=data,
                bounce_quantization=uncertainty,local_symplectic_audits=audits)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", "--study", dest="study", type=Path,required=True, help="Prior study root with hilltop/ and background/")
    parser.add_argument("--archive",type=Path, help="Default is STUDY/archive")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--detunings",type=float,nargs="+",default=[1e-3,1e-10,1e-20])
    parser.add_argument("--nmax",type=int,default=16)
    parser.add_argument("--rtol",type=float,default=1e-11)
    parser.add_argument("--max-step",type=float,default=.035)
    args = parser.parse_args()
    if args.nmax < 2:
        parser.error("--nmax must be at least two")
    M = module(args.study/"hilltop"/"hilltop.py", "parent_hill")
    C = module(args.study/"hilltop"/"check.py", "parent_hill_check")
    mod, ref = M.BG.source(args.archive or args.study/"archive")
    tables = mod.build_tables()
    history = json.loads((args.study/"hilltop"/"HILLTOP.json").read_text())
    branches = {}
    for key,branch in history["branches"].items():
        p = M.BG.Potential(mod,tables,float(key.split("=")[1]),ref[key])
        hill = M.ClosedLaunch(p)
        force,_ = C.local_force(p,hill.top)
        ds = {}
        for det,record in branch["detunings"].items():
            if float(det) in args.detunings:
                ds[det]=integrate(p,hill,force,record["primary"]["log_abs_bounce_seed"],
                                 np.arange(2,args.nmax+1),[1e-5,1e-4,1e-3,1e-2],args.rtol,args.max_step,
                                 float(det),record["primary"]["x_match"])
                print(key,det, [row["max_delta_theta_row_change"] for row in ds[det]["endpoints"]], flush=True)
        branches[key] = ds
    args.output.write_text(json.dumps(dict(branches=branches),indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
