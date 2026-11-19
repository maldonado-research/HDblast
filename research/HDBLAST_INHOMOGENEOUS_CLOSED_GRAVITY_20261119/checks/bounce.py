"""Short-time proper-mean derivatives evaluated directly on the initial slice."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval
from scipy.special import roots_legendre


def load_branches(previous):
    path = Path(__file__).resolve().parents[1] / "initial/constraint.py"
    spec = importlib.util.spec_from_file_location("initial_constraint_for_checks", path)
    initial = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(initial)
    return initial.load_branches(previous)


def force_derivative(branch, q):
    """Differentiate the actual quotient on each spline segment."""
    p = branch.fitted
    field = branch.top + q
    phi, f = p.spline(field).T
    dphi, df = p.spline(field, 1).T
    ddphi, ddf = p.spline(field, 2).T
    z = 1 - phi / p.pm
    weight = 1 + p.b * phi
    n = z**2 * weight
    n_phi = -2 * z * weight / p.pm + z**2 * p.b
    n_phi_phi = 2 * weight / p.pm**2 - 4 * z * p.b / p.pm
    n1 = n_phi * dphi
    n2 = n_phi_phi * dphi**2 + n_phi * ddphi
    return -3 * (n2 - 4 * n1 * df / f + n * (6*(df/f)**2 - 2*ddf/f)) / (f**2*p.v0)


def derivatives(branch, amplitude, offset, coefficients, count):
    # A Gauss–Legendre rule in chi is independent of the Galerkin projection grid.
    nodes, weights = roots_legendre(count)
    chi = np.pi * (nodes+1)/2
    z = np.cos(chi)
    sin2 = np.sin(chi)**2
    # U_(2j) = 1 + 2(T_2 + ... + T_(2j)).
    cosine = np.zeros(2*len(coefficients)-1)
    cosine[0] = np.sum(coefficients)
    cosine[2::2] = 2*np.cumsum(coefficients[:0:-1])[::-1]
    psi = chebval(z, cosine)
    psi_z = chebval(z, chebder(cosine))
    psi_zz = chebval(z, chebder(cosine, 2))
    sigma_z = 2*psi_z/psi
    sigma_zz = 2*psi_zz/psi - sigma_z**2/2
    scale = branch.h**2 / psi**4
    proper_weight = weights*sin2*psi**6
    proper_weight /= sum(proper_weight)
    avg = lambda f: float(proper_weight @ f)
    cov = lambda f, g: avg((f-avg(f))*(g-avg(g)))

    quadrupole = amplitude*(4*z**2-1)
    q = offset + quadrupole
    q_z = 8*amplitude*z
    gradient2 = scale*sin2*q_z**2
    lap = scale*(-8*quadrupole + sin2*sigma_z*q_z)
    u, force = branch.potential(q)
    q2 = lap + force
    # Ricci eigenvalues of gamma = exp(2 sigma) * round_metric.
    ric_r = scale*(2-2*sin2*sigma_zz+4*z*sigma_z)
    ric_t = scale*(2-sin2*sigma_zz+4*z*sigma_z-sin2*sigma_z**2)
    hr1 = 3*u + gradient2 - ric_r
    ht1 = 3*u - ric_t
    theta1 = 3*u
    theta3 = -force*q2 - 2*q2**2 - 2*(hr1**2 + 2*ht1**2)
    jac4 = theta3 + 3*theta1**2
    m0 = avg(q)
    m2 = avg(force) + cov(q, theta1)
    m4 = (avg(force_derivative(branch,q)*q2) + avg(theta1*(4*q2-lap))
          + cov(q, jac4) - 6*avg(theta1)*m2)
    v2 = avg(theta1)/3
    v4 = avg(jac4)/3 - 2*avg(theta1)**2/3
    return dict(mean=m0, mean_xx=m2, mean_xxxx=m4,
                relative_volume_radius_xx=v2, relative_volume_radius_xxxx=v4,
                max_initial_hamiltonian=float(np.max(abs(ric_r+2*ric_t-gradient2-6*u))),
                mean_laplacian=avg(lap),
                early_times={str(x):dict(mean_quadratic=m0+m2*x*x/2,
                                         mean_quartic=m0+m2*x*x/2+m4*x**4/24)
                             for x in (.01,.02)})


def compare_to_evolution(prediction, run):
    errors = []
    for name, early in prediction["early_times"].items():
        actual = run["samples"][name]["mean"]
        early["evolved_mean"] = actual
        early["quadratic_error"] = early["mean_quadratic"]-actual
        error = early["mean_quartic"]-actual
        early["quartic_error"] = error
        change = run["samples"][name]["mean_change"]
        early["quartic_error_relative_to_change"] = error/change if change else None
        errors.append(error)
    ratio = errors[1]/errors[0] if errors[0] else None
    if run["amplitude"] and (ratio is None or abs(ratio-64)>1):
        raise ArithmeticError("Fourth-order Taylor error fails its early-time scaling check")
    prediction["doubled_time_quartic_error_ratio"] = ratio


def calculate(previous, initial_file, evolution_file):
    data = json.loads(initial_file.read_text())
    evolved = json.loads(evolution_file.read_text())["branches"]
    output = {}
    for key, branch in load_branches(previous).items():
        output[key] = {}
        for group in ("runs", "proper_mean_centered_runs"):
            rows = {}
            for name, run in evolved[key][group].items():
                amp = run["amplitude"]
                record = data["branches"][key][group][str(float(amp))]
                coefficients = record["psi_even_U_coefficients"]
                args = (branch, amp, run["offset"], coefficients)
                a = derivatives(*args, 512)
                b = derivatives(*args, 768)
                for field in ("mean", "mean_xx", "mean_xxxx"):
                    if abs(a[field]-b[field]) > 3e-11:
                        raise ArithmeticError(f"Quadrature check failed for {key} {name} {field}")
                if abs(b["mean_xx"]-record["proper_mean_q_acceleration"]) > 2e-12:
                    raise ArithmeticError("Proper-mean acceleration disagrees with the initial constraint results")
                if b["max_initial_hamiltonian"] > 1e-11 or abs(b["mean_laplacian"]) > 1e-14:
                    raise ArithmeticError("Independent initial geometry check failed")
                compare_to_evolution(b, run)
                if group == "runs" and amp:
                    refinement = next(item for item in evolved[key]["spatial_refinements"]
                                      if item["maximum_degree"]==24 and item["amplitude"]==amp)
                    errors = [entry["mean_quartic"]-refinement["samples"][time]["mean"]
                              for time, entry in b["early_times"].items()]
                    b["degree_24_doubled_time_quartic_error_ratio"] = errors[1]/errors[0]
                    if abs(errors[1]/errors[0]-64) > 1:
                        raise ArithmeticError("Degree 24 evolution fails its Taylor check")
                rows[name] = b
            output[key][group] = rows
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--initial", type=Path, required=True)
    parser.add_argument("--evolution", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    result = calculate(options.previous, options.initial, options.evolution)
    options.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
