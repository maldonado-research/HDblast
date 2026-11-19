"""Guarded continuation of the pinned closed-sphere Einstein–scalar evolutions."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
from scipy.integrate import DOP853
from scipy.optimize import brentq, minimize_scalar
from scipy.special import eval_chebyu

FIELD_MARGIN = 1e-4
RESIDUAL_LIMIT = 1e-4
AMPLITUDES = (.001, -.001, .0005, -.0005)
CONFIG = {"coarse": (64, .003), "medium": (112, .003), "fine": (192, .003),
          "time": (192, .0015), "higher": (256, .003), "near_edge": (192, .003)}


class DomainStep(Exception):
    """A trial RK stage lies outside the fitted domain."""


def load(gravity):
    spec = importlib.util.spec_from_file_location("pinned_evolution", gravity/"evolution/evolve.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Monitor:
    def __init__(self, solver, amplitude, margin):
        self.margin = margin
        self.e = solver
        self.b = solver.branch
        self.sign = self.b.fitted.sign
        self.target = self.b.fitted.tmin-self.b.top
        self.low, self.high = self.b.fitted.theta[[0, -1]]-self.b.top
        # Reflection symmetry identifies the other hemisphere exactly.
        mid = (solver.monitor_chi.size+1)//2
        self.chi = solver.monitor_chi[:mid]
        self.y = np.cos(self.chi)
        self.w = 1-self.y**2
        self.mat = tuple(mat[:mid] for mat in solver.monitor_all)
        self.U, self.Uy, self.Uyy = self.mat
        # symmetry-fixed location of maximum displacement toward the minimum at launch
        n = solver.degrees
        self.curvature_weights = n*(n+1)*(n+2)/3
        self.i = 0 if self.sign*amplitude>0 else -1

    def certify_domain(self, c, values):
        """Enclose the full angular polynomial, including gaps between samples.

        U_n(cos chi) is a sum of Fourier modes -n,-n+2,...,n. Their
        squared frequencies sum to n(n+1)(n+2)/3. The resulting bound on
        |q''(chi)| bounds linear-interpolation error by B*spacing^2/8.
        """
        bound = (abs(c)@self.curvature_weights)*(self.chi[1]-self.chi[0])**2/8
        if values.min()-bound>self.low and values.max()+bound<self.high:
            return
        left, right = self.chi[:-1], self.chi[1:]
        ql, qr = values[:-1], values[1:]
        for _ in range(16):
            if min(ql.min(), qr.min())<=self.low or max(ql.max(), qr.max())>=self.high:
                raise DomainStep
            uncertain = ((np.minimum(ql, qr)-bound<=self.low) |
                         (np.maximum(ql, qr)+bound>=self.high))
            if not uncertain.any():
                return
            left, right, ql, qr = (a[uncertain] for a in (left, right, ql, qr))
            mid = (left+right)/2
            qm = eval_chebyu(self.e.degrees, np.cos(mid)[:, None])@c
            left, right = np.r_[left, mid], np.r_[mid, right]
            ql, qr = np.r_[ql, qm], np.r_[qm, qr]
            bound /= 4
        raise DomainStep

    def guarded_equation(self, t, z):
        c = z.reshape(6, -1)[4]
        qm = self.U@c
        self.certify_domain(c, qm)
        qc = self.e.U@c
        if min(qc.min(), qm.min()) <= self.low or max(qc.max(), qm.max()) >= self.high:
            raise DomainStep
        return self.e.equation(t, z)

    def inspect(self, t, z):
        e = self.e
        fields = e.fields(z, self.mat)
        (m, d, b, shear, q, pi), (my, dy, by, qy), (myy, qyy) = fields
        a, H, inv, ly, F, dK, dJ, dric, dR = e.geometry(t, fields, self.y, self.w)
        u, force = self.b.potential(q)
        beta = H+b
        hlocal = beta+self.w*shear/3
        ham = (dR+12*H*b+6*b*b+4*self.w*beta*shear-pi*pi
               -inv*self.w*qy*qy-6*(u-self.b.u_top))
        mom = np.sqrt(inv*self.w)*(2*(by+F*shear)+pi*qy)
        hc = float(max(abs(ham))/(6*self.b.u_top))
        mc = float(max(abs(mom))/self.b.u_top)
        # Independent higher-order Gauss-U quadrature for proper-volume observables
        M, D, B, E = (e.check[0] @ z.reshape(6, -1)[:4].T).T
        measure = e.check_weight*np.exp(3*M+e.check_w*D)
        Hv = H+measure@(B+e.check_w*E/3)/measure.sum()
        volume = 2*np.pi**2*a**3*measure.sum()
        imin = int(np.argmin(hlocal))
        if self.sign>0:
            iq = int(np.argmax(q))
        else:
            iq = int(np.argmin(q))
        i = self.i
        lap = inv[i]*(self.w[i]*qyy[i]-3*self.y[i]*qy[i]+self.w[i]*(2*my[i]-ly[i])*qy[i])
        damping = -3*hlocal[i]*pi[i]
        s = self.sign
        # Use scalars oriented so each physical event changes positive to negative.
        return dict(q=q, pi=pi, hlocal=hlocal, iq=iq, ih=imin,
                    monitored=dict(domain_guard=min(q.min()-self.low, self.high-q.max())-self.margin,
                                   constraint_guard=RESIDUAL_LIMIT-max(hc, mc),
                                   first_minimum=s*(self.target-q[iq]),
                                   first_local_contraction=hlocal[imin], volume_turnaround=Hv,
                                   return_through_minimum=s*(q[self.i]-self.target),
                                   first_turn_after_entry=s*pi[self.i],
                                   second_turn_after_entry=-s*pi[self.i],
                                   third_through_minimum=s*(self.target-q[self.i])),
                    coeff=z.reshape(6,-1)[4], hc=hc, mc=mc, Hv=Hv, volume=volume, curvature_bound=float(abs(z.reshape(6,-1)[4])@self.curvature_weights),
                    lap=float(lap), damping=float(damping), force=float(force[i]), Jmin=float(min(np.exp(3*m+self.w*d))))

    def snapshot(self, t, data):
        q, pi, h = data["q"], data["pi"], data["hlocal"]
        iq, ih = data["iq"], data["ih"]
        i = self.i
        def extremum(sign):
            index = int(np.argmax(sign*q))
            if index in (0,len(q)-1):
                return float(q[index])
            def negative_value(chi):
                return -sign*(eval_chebyu(self.e.degrees, np.cos(chi))@data["coeff"])
            best = minimize_scalar(negative_value, bounds=(self.chi[index-1],self.chi[index+1]),
                                   method="bounded", options={"xatol":1e-13})
            return float(-sign*best.fun)
        return dict(x=float(t), q_min=extremum(-1), q_max=extremum(1),
                    chi_extreme_in_roll_direction=float(self.chi[iq]), H_at_extreme_in_roll_direction=float(h[iq]),
                    H_local_min=float(h[ih]), chi_H_local_min=float(self.chi[ih]),
                    H_volume=float(data["Hv"]), volume=float(data["volume"]),
                    min_spatial_volume_over_reference=data["Jmin"],
                    ham_max=data["hc"], momentum_max=data["mc"], first_worldline_q=float(q[i]),
                    first_worldline_Pi=float(pi[i]), first_worldline_H=float(h[i]),
                    first_worldline_acceleration_terms={"laplacian":data["lap"],
                                          "expansion":data["damping"], "potential_force":data["force"]},
                    field_extremum_interpolation_error_bound=float(data["curvature_bound"]*
                                                        (self.chi[1]-self.chi[0])**2/8))


def run(e, amplitude, old, dt, rtol, margin):
    m = Monitor(e, amplitude, margin)
    z = e.initial(amplitude, 0.)
    t = 0.
    dt_initial = dt
    sol = None
    events = {}
    maxima = np.zeros(2)
    last = m.inspect(t, z)
    n = 0
    retries = 0
    crossed = np.zeros(len(m.chi), dtype=int)
    qold = last["q"]
    old_x = old["x"]
    cap_delta = None
    # Bounds on every potential evaluation are checked at the RHS. If a trial
    # stage hits the hard boundary, discard it and restart from the accepted step.
    while t < 1.1:
        if sol is None:
            sol = DOP853(m.guarded_equation, t, z, 1.1, max_step=dt,
                         rtol=rtol, atol=1e-14 if rtol>1e-11 else 2e-15)
        try:
            sol.step()
            if sol.status == "failed":
                raise ArithmeticError("DOP853 step failed")
            dense = sol.dense_output()
        except DomainStep:
            smaller = min(dt, sol.h_abs, sol.t-t if sol.t>t else dt)/2
            if smaller < 1e-10:
                raise ArithmeticError("Insufficient spline guard margin")
            sol = DOP853(m.guarded_equation, t, z, 1.1, first_step=smaller,
                         max_step=smaller, rtol=rtol, atol=1e-14 if rtol>1e-11 else 2e-15)
            dt = smaller
            retries += 1
            continue
        end = sol.t
        new = m.inspect(end, sol.y)
        def at(x):
            return m.inspect(x, dense(x))
        def find(name):
            return brentq(lambda x: at(x)["monitored"][name], t, end, xtol=5e-14)
        candidates = {}
        for name in last["monitored"]:
            if name in events:
                continue
            forerunner = {"first_turn_after_entry":"first_minimum",
                          "return_through_minimum":"first_turn_after_entry",
                          "second_turn_after_entry":"return_through_minimum",
                          "third_through_minimum":"second_turn_after_entry"}.get(name)
            if forerunner and forerunner not in events:
                continue
            if last["monitored"][name]>0 and new["monitored"][name]<=0:
                candidates[name] = find(name)
        terminal = min([candidates[k] for k in ("domain_guard", "constraint_guard")
                        if k in candidates], default=None)
        for name, x in candidates.items():
            if terminal is None or x<=terminal:
                events[name] = m.snapshot(x, at(x))
        if terminal is not None:
            end = terminal
            new = at(end)
        if t < old_x <= end:
            zp = dense(old_x)
            r = e.record(old_x, zp)
            cap_delta = {k:float(r[k]-old[k]) for k in
                         ("maximum_field", "mean", "H_volume", "volume")}
        crossed += ((qold-m.target)*(new["q"]-m.target) < 0)
        qold = new["q"]
        maxima = np.maximum(maxima, [new["hc"], new["mc"]])
        n += 1
        t, z, last = end, sol.y, new
        if terminal is not None:
            break
    stop = next((k for k in ("constraint_guard","domain_guard") if k in events), "time_cap")
    return dict(degree=e.degree, nodes=len(e.chi), field_guard_margin=margin, max_step=dt_initial, rtol=rtol, monitor_hemisphere_nodes=len(m.chi),
                time_steps=n, hard_domain_retries=retries, first_worldline_chi=float(m.chi[m.i]),
                stop=stop, events=events, end=m.snapshot(t, last),
                max_crossings_on_any_sampled_worldline=int(max(crossed)),
                old_cap_x=old_x, old_cap_differences=cap_delta,
                max_step_sampled_ham=float(maxima[0]), max_step_sampled_momentum=float(maxima[1]))


def calculate(gravity, levels, amps):
    module = load(gravity)
    old = json.loads((gravity/"evolution/results.json").read_text())["branches"]
    result = {"field_guard_margin":FIELD_MARGIN,"constraint_limit":RESIDUAL_LIMIT,"branches":{}}
    for key, branch in module.constraint.load_branches(gravity/"previous").items():
        entry = result["branches"][key] = dict(
            canonical_minimum_offset=branch.fitted.tmin-branch.top,
            spline_offset_bounds=(branch.fitted.theta[[0,-1]]-branch.top).tolist(),
            residual_slope=branch.du_top.item(), runs={})
        for level in levels:
            degree, dt = CONFIG[level]
            e = module.Evolution(branch,degree,2*degree+32,monitor_nodes=6*degree+1)
            for A in amps:
                a = f"{A:g}"
                old_stop = old[key]["controlled_extensions"]["runs"][a]["at_stop"]
                r = run(e, A, old_stop, dt, 3e-12 if level=="time" else 3e-11,
                        1e-6 if level=="near_edge" else FIELD_MARGIN)
                entry["runs"].setdefault(a, {})[level] = r
                print(key, a, level, r["stop"], round(r["end"]["x"],7),
                      {k:round(v["x"],7) for k,v in r["events"].items()},flush=True)
    return result


def summarize(data):
    """Verify scientific events and retain small, comparable refinement records."""
    max_time = max_spatial = max_old = max_edge = max_h_delta = max_q_delta = 0.
    for key, entry in data["branches"].items():
        for label, levels in entry["runs"].items():
            if set(levels) != set(CONFIG):
                continue
            a = float(label)
            primary = levels["fine"]
            names = set(primary["events"])
            if primary["stop"]!="domain_guard" or "first_local_contraction" not in names:
                raise ArithmeticError("Incomplete physical event record")
            neg = key=="phi_m=-0.5"
            if ("first_minimum" in names)!=neg:
                raise ArithmeticError("Wrong first-passage classification")
            oscillating = neg and a>0
            if ("return_through_minimum" in names)!=oscillating:
                raise ArithmeticError("Wrong return-passage classification")
            turnaround = oscillating or (not neg and a==-.001)
            if ("volume_turnaround" in names)!=turnaround:
                raise ArithmeticError("Wrong volume-turnaround classification")
            count = 2 if oscillating else int(neg)
            if primary["max_crossings_on_any_sampled_worldline"]!=count:
                raise ArithmeticError("Inconsistent worldline crossing count")
            if neg:
                first = primary["events"]["first_minimum"]
                if (abs(first["chi_extreme_in_roll_direction"]-primary["first_worldline_chi"])>1e-10
                        or min(first["first_worldline_H"],first["H_volume"],
                               first["first_worldline_Pi"])<=0):
                    raise ArithmeticError("First passage location or direction failed")
            if oscillating:
                at_turn = primary["events"]["first_turn_after_entry"]
                at_return = primary["events"]["return_through_minimum"]
                if at_return["first_worldline_Pi"]>=0 or abs(at_turn["first_worldline_Pi"])>1e-7:
                    raise ArithmeticError("Unresolved outward turn and return")
            for name in ("time","higher","near_edge"):
                r = levels[name]
                if set(r["events"])!=names or r["max_crossings_on_any_sampled_worldline"]!=count:
                    raise ArithmeticError("Event sequence changed under refinement")
                for key_ev in names:
                    d = abs(r["events"][key_ev]["x"]-primary["events"][key_ev]["x"])
                    if name=="time":
                        max_time=max(max_time,d)
                    elif name=="higher":
                        max_spatial=max(max_spatial,d)
                    else:
                        max_edge=max(max_edge,d)
                if name=="higher":
                    max_h_delta = max(max_h_delta, abs(r["end"]["H_volume"]-primary["end"]["H_volume"]))
                    max_q_delta = max(max_q_delta, *(abs(r["end"][f]-primary["end"][f])
                                                      for f in ("q_min","q_max","first_worldline_q")))
            max_old = max(max_old, *(abs(v) for v in primary["old_cap_differences"].values()))
            for name, r in levels.items():
                if name=="fine":
                    continue
                r["event_times"] = {event: v["x"] for event,v in r.pop("events").items()}
                r["end"] = {k:r["end"][k] for k in ("x","q_min","q_max","H_volume",
                          "first_worldline_q","first_worldline_Pi","ham_max","momentum_max")}
    if max(max_time,max_spatial,max_q_delta,max_h_delta)>3e-7 or max_old>2e-9 or max_edge>2e-5:
        raise ArithmeticError(f"Failed convergence or regression: time={max_time:g}, spatial={max_spatial:g}, field={max_q_delta:g}, H={max_h_delta:g}, old={max_old:g}, guard={max_edge:g}")
    if all(set(r)==set(CONFIG) for b in data["branches"].values() for r in b["runs"].values()):
        data["verification"] = dict(max_timestep_event_time_difference=max_time,
                       max_higher_degree_event_time_difference=max_spatial,
                       max_higher_degree_endpoint_H_volume_difference=max_h_delta,
                       max_higher_degree_endpoint_field_difference=max_q_delta,
                       max_prior_cap_observable_difference=max_old,max_guard_event_time_difference=max_edge)
    return data


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--gravity",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--levels",nargs="+",choices=CONFIG,default=list(CONFIG))
    p.add_argument("--amplitudes",type=float,nargs="+",default=AMPLITUDES)
    args = p.parse_args()
    result = summarize(calculate(args.gravity, args.levels,args.amplitudes))
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
