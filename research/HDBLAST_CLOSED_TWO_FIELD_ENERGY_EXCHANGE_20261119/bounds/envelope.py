"""Continuum energy inequalities evaluated on the archived first-passage spacetimes."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, simpson
from scipy.optimize import minimize_scalar
from scipy.special import eval_chebyu

HERE = Path(__file__).resolve().parent


def module_at(path):
    spec = importlib.util.spec_from_file_location("canonical_passage_for_envelope", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(geometry, solution, cut, angular_points, time_points):
    # Equatorial reflection symmetry reduces the search to half the sphere.
    angle = np.linspace(0., np.pi/2, angular_points)
    degrees = geometry.degrees
    basis = eval_chebyu(degrees, np.cos(angle)[:, None])
    t = np.linspace(cut-.005, cut+.005, time_points)
    # min/max of local instantaneous growth and |d(log Z)/dt|
    rates = []
    unpolished_error = 0.
    for time in t:
        state = solution.sol(time).reshape(6, geometry.n)
        H = geometry.background(time)[1]
        def coefficients(chi, U):
            b, e, q, pi = (U @ state[[2, 3, 4, 5]].T).T
            Ht = H+b
            Hr = Ht+np.sin(chi)**2*e
            expansion = Hr+2*Ht
            spline = geometry.branch.fitted.spline
            field = geometry.branch.top+q
            s = -pi*spline(field, 1)[..., 1]/spline(field)[..., 1]
            return np.asarray([-expansion-s, expansion-2*Hr+s, expansion-2*Ht+s, s])
        columns = coefficients(angle, basis)
        def extremum(row, direction):
            values = direction*columns[row]
            candidates = [float(values[0]), float(values[-1])]
            for i in np.where((values[1:-1]>values[:-2]) & (values[1:-1]>=values[2:]))[0]+1:
                def negative(theta):
                    U = eval_chebyu(degrees, np.cos(theta))
                    return float(-direction*coefficients(theta,U)[row])
                answer = minimize_scalar(negative, bounds=(angle[i-1], angle[i+1]),
                                         method="bounded", options={"xatol":1e-11})
                candidates.append(-answer.fun)
            return max(max(candidates), float(np.max(values)))
        lower = -max(extremum(i, -1) for i in range(3))
        upper = max(extremum(i, 1) for i in range(3))
        absolute_s = max(extremum(3,-1), extremum(3,1))
        unpolished_error = max(unpolished_error, np.min(columns[:3])-lower,
                               upper-np.max(columns[:3]), absolute_s-np.max(abs(columns[3])))
        rates.append([lower, upper, absolute_s])

    rates = np.asarray(rates)
    int_low, int_high = [cumulative_simpson(column,x=t,initial=0.) for column in rates.T[:2]]
    cap = simpson(rates[:,2]*np.exp(int_high), x=t)
    return dict(begin=float(t[0]), end=float(t[-1]), E_end_over_E_begin_lower=float(np.exp(int_low[-1])),
                E_end_over_E_begin_upper=float(np.exp(int_high[-1])),
                abs_coupling_work_over_E_begin_upper=float(cap),
                integrated_max_abs_dlogZ=float(simpson(rates[:,2], x=t)),
                largest_spatial_polish=float(unpolished_error), center_rates=rates[len(t)//2].tolist())


def calculate(prior):
    module = module_at(prior/"passage/passage.py")
    solver = module.dynamic.intrinsic.gravity_module(prior/"gravity")
    branch = solver.constraint.load_branches(prior/"gravity/previous")["phi_m=-0.5"]
    source = json.loads((prior/"passage/PASSAGE.json").read_bytes())["runs"]
    reports = {}
    for amp in (.001,-.001):
        rows = {}
        for label, degree, nodes, angular, times in [("coarse",48,192,121,33),("fine",80,256,241,65)]:
            g, trajectory, info = module.geometry_through_passage(solver,branch,amp,degree,nodes)
            if abs(info["crossing_x"]-source[f"{amp:g}"]["crossing_x"]) > 2e-9:
                raise ArithmeticError("Archived event time mismatch")
            rows[label] = evaluate(g, trajectory, info["crossing_x"], angular,times)
        errors = [abs(rows['coarse'][k]-rows['fine'][k]) for k in rows['fine']
                  if k not in ('center_rates','largest_spatial_polish')]
        rows['largest_bound_refinement'] = max(errors)
        if max(errors) >= 3e-7:
            raise ArithmeticError(f"Envelope refinement error {max(errors)}")
        reports[f"{amp:g}"] = rows
    return reports


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, default=HERE.parent/"prior")
    parser.add_argument("--output", type=Path, required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior.resolve()),indent=2,allow_nan=False)+"\n")
