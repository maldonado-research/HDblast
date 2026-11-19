"""Closed de Sitter scalar transfer and finite Gaussian band diagnostics."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import hyp2f1, log_ndtr

SOURCE_SHA = "f89850b46c6a93303bec3968c4712659796f53fc8ac20fad593342c6b5db86b6"
VOLUME = 2 * math.pi**2
SHELLS = np.arange(17)
K2 = SHELLS * (SHELLS + 2)
DEGENERACY = (SHELLS + 1)**2
FILTERS = {
    "all_2_to_16": (SHELLS >= 2).astype(float),
    "band_2_to_8": ((SHELLS >= 2) & (SHELLS <= 8)).astype(float),
    "heat_angle_0.25": (SHELLS >= 2) * np.exp(-0.5 * 0.25**2 * K2),
    "heat_angle_0.5": (SHELLS >= 2) * np.exp(-0.5 * 0.5**2 * K2),
    "quadrupole_n2": (SHELLS == 2).astype(float),
    "dipole_n1_control": (SHELLS == 1).astype(float),
}


def prior_model(prior):
    spec = importlib.util.spec_from_file_location("preparation_hilltop", prior / "hilltop/hilltop.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    model, ref = module.BG.source(prior / "archive")
    return module, model, ref, model.build_tables()


def force_radii(p, top, mu2):
    """First 10% crossing on each side within the tabulated neighborhood."""
    radii = []
    for sign in [-1, 1]:
        def residual(d):
            linear = mu2 * sign * d
            return abs((-3 * p.at(top + sign*d)[1] - linear) / linear) - 0.10
        grid = np.linspace(1e-5, 0.25, 1001)
        lo = grid[0]
        for hi in grid[1:]:
            if residual(hi) > 0:
                radii.append(brentq(residual, lo, hi))
                break
            lo = hi
        else:
            raise ValueError("10% hilltop force radius missing")
    return radii


def transfer(x, h, mu2):
    """[n, output q/v, input q/v], analytic even and odd solutions at all n."""
    t = h*x
    u = math.tanh(t)
    sech2 = math.cosh(t)**-2
    disc = math.sqrt(9 + 4*mu2/h**2)
    growth = 0.5 * (disc - 3)
    common = math.cosh(t)**growth
    # Pfaff transformation keeps the hypergeometric argument on [0,1].
    a = SHELLS/2 + (3 + disc)/4
    b = SHELLS/2 + (3 - disc)/4
    matrix = np.empty((len(SHELLS), 2, 2))
    for j in [0, 1]:
        A, B, C = 0.5 - a + j/2, b + j/2, 0.5 + j
        F = hyp2f1(A, B, C, u*u)
        Fz = A*B/C * hyp2f1(A+1, B+1, C+1, u*u)
        prefactor = common * (u/h)**j
        prefactor_prime = growth*h*u*prefactor + j*common*sech2
        matrix[:, 0, j] = prefactor * F
        matrix[:, 1, j] = prefactor_prime*F + prefactor * 2*h*u*sech2*Fz
    return matrix


def integrate(start, end, h, mu2, dense):
    def equation(x, flat):
        f = flat.reshape(len(SHELLS), 2, 2)
        d = np.empty_like(f)
        d[:, 0] = f[:, 1]
        d[:, 1] = ((mu2 - K2*h*h/math.cosh(h*x)**2)[:, None]*f[:, 0]
                   - 3*h*math.tanh(h*x)*f[:, 1])
        return d.ravel()
    eye = np.tile(np.eye(2), (len(SHELLS), 1, 1))
    result = solve_ivp(equation, (start, end), eye.ravel(), method="DOP853",
                       rtol=1e-12, atol=1e-14, max_step=0.04, dense_output=dense)
    if not result.success:
        raise ArithmeticError(result.message)
    return result


def validation(times, h, mu2):
    sol = integrate(0, max(times), h, mu2, True)
    transfer_error = 0.
    wronskian_error = 0.
    for x in times:
        exact = transfer(x, h, mu2)
        ode = sol.sol(x).reshape(len(SHELLS), 2, 2)
        transfer_error = max(transfer_error, np.max(abs(exact - ode) / np.maximum(1., abs(exact))))
        # Compare local determinants across the entire domain without subtracting
        # products of nearly parallel exponentially growing columns at late times.
        x1 = x + 0.07
        local = integrate(x, x1, h, mu2, False).y[:, -1].reshape(len(SHELLS), 2, 2)
        abel = (math.cosh(h*x) / math.cosh(h*x1))**3
        wronskian_error = max(wronskian_error, np.max(abs(np.linalg.det(local)/abel - 1)))
    early = np.linspace(0., 0.3, 16)
    early_error = max(np.max(abs(np.linalg.det(transfer(x, h, mu2))*math.cosh(h*x)**3 - 1))
                      for x in early)
    if transfer_error > 5e-8 or max(early_error, wronskian_error) > 1e-9:
        raise ArithmeticError((transfer_error, early_error, wronskian_error))
    return dict(max_entry_error_scaled_by_max_1_abs_analytic=transfer_error,
                max_early_analytic_relative_Wronskian_error=early_error,
                max_local_ODE_relative_Wronskian_error=wronskian_error,
                local_Wronskian_interval=0.07, analytic_early_end=0.3)


def moments(x, h, mu2, H_over_M):
    f = transfer(x, h, mu2)
    omega = np.sqrt(1 + h*h*K2)
    # Each orthonormal real mode is prepared with no field-velocity covariance.
    # Summing over degeneracies gives pointwise moments of each independent shell.
    normalization = DEGENERACY/VOLUME * (H_over_M**2 * h**3/2)
    qq = normalization * (f[:, 0, 0]**2 / omega + omega*f[:, 0, 1]**2)
    vv = normalization * (f[:, 1, 0]**2 / omega + omega*f[:, 1, 1]**2)
    qv = normalization * (f[:, 0, 0]*f[:, 1, 0] / omega + omega*f[:, 0, 1]*f[:, 1, 1])
    grad = qq * K2 * h*h / math.cosh(h*x)**2
    return qq, vv, qv, grad


def first_crossing(function, limit, grid, values):
    if values[0] >= limit:
        return 0.
    for i in range(1, len(grid)):
        if values[i] >= limit:
            return brentq(lambda x: function(x)-limit, grid[i-1], grid[i], xtol=1e-13)
    raise ArithmeticError("Crossing missing")


def case_result(h, mu2, data, radii):
    H_over_M = 2*math.pi*data["archive_displacement_scale"]
    primary = data["primary"]
    residence = primary["x_match"]
    # Dense early sampling detects possible oscillatory threshold crossings.
    grid = np.linspace(0., residence, math.ceil(residence / .008) + 1)
    samples = np.array([moments(x, h, mu2, H_over_M) for x in grid])
    out = dict(H_top_over_Mpl=H_over_M, homogeneous_linear_handoff_x=residence,
               homogeneous_linear_handoff_abs_delta=1e-5, filters={})
    kappa_field = H_over_M**2 * h**3 / VOLUME
    def minimax_residual(x):
        return 9 * kappa_field * transfer(x, h, mu2)[2, 0, 1] / 2 - min(radii)**2
    out["x_n2_optimal_spatial_rms_at_force_radius"] = brentq(minimax_residual, 0., residence)
    for label, window in FILTERS.items():
        weights = window**2
        def field(x):
            return math.sqrt(moments(x, h, mu2, H_over_M)[0] @ weights)
        rms = np.sqrt(samples[:, 0, :] @ weights)
        dx = first_crossing(field, min(radii), grid, rms)
        overlap = first_crossing(field, 1e-5, grid, rms)
        out["filters"][label] = dict(bounce_variance=rms[0]**2,
                                      x_rms_equals_1e_minus5=overlap,
                                      x_rms_at_force_radius=dx,
                                      homogeneous_handoff_x_minus_force_crossing=residence-dx)
    weights = FILTERS["all_2_to_16"]
    def energies(q, v, grad):
        # Densities divided by stationary V*=3 Mpl² Htop² h².
        return np.stack((v, grad, -mu2*q),axis=-1) / (6*h*h)
    energy_series = energies(samples[:, 0, :] @ weights, samples[:, 1, :] @ weights,
                             samples[:, 3, :] @ weights)
    def energy_components(x):
        qq, vv, _, grad = moments(x, h, mu2, H_over_M)
        return energies(qq@weights, vv@weights, grad@weights)
    def largest_energy(x):
        return np.max(abs(energy_components(x)))
    x_energy = first_crossing(largest_energy, .1, grid, np.max(abs(energy_series), axis=1))
    x_force = out["filters"]["all_2_to_16"]["x_rms_at_force_radius"]
    for name, x in [("bounce", 0.), ("at_force_crossing", x_force), ("at_energy_crossing", x_energy)]:
        K, G, U = energy_components(x)
        out[name+"_quadratic_energy_over_Vstar"] = dict(kinetic=K, gradient=G, potential=U,
                                                        signed_sum=K+G+U)
    out["x_largest_abs_energy_component_equals_0.1_Vstar"] = x_energy
    x_diagnostic = min(x_force, x_energy)
    out["first_of_two_full_band_diagnostics_x"] = x_diagnostic
    out["variance_snapshots"] = []
    for x in [0., 0.5*x_diagnostic, x_diagnostic]:
        q, _, _, _ = moments(x, h, mu2, H_over_M)
        variances = {label: float(q @ (window**2)) for label, window in FILTERS.items()}
        sigma = math.sqrt(variances["all_2_to_16"])
        log_tail = np.logaddexp(*(log_ndtr(-d/sigma) for d in radii))
        out["variance_snapshots"].append(dict(
            x=x, point_variance_by_filter=variances,
            raw_band_variance_fraction_by_n_2_to_16=(q[2:]/sum(q[2:])).tolist(),
            log10_one_point_gaussian_weight_outside_asymmetric_force_window=log_tail/math.log(10)))
    return out


def run(prior, destination):
    raw = (prior/"hilltop/HILLTOP.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA:
        raise ValueError("Unexpected HILLTOP.json")
    history = json.loads(raw)
    module, model, ref, tables = prior_model(prior)
    result = dict(input_hilltop_sha256=SOURCE_SHA, maximum_prepared_n=16,
                  S3_unit_volume=VOLUME, reference_mass_over_Htop=1., branches={})
    matrices = []
    for key, branch in history["branches"].items():
        p = module.BG.Potential(model, tables, float(key.split("=")[1]), ref[key])
        top, h, mu2, _ = module.field_top(p)
        if abs(mu2 - branch["negative_mass_squared_over_Htop_squared"]) > 1e-8:
            raise ArithmeticError("Archived mass mismatch")
        radii = force_radii(p, top, mu2)
        times = sorted({0., .1, .25, .5, 1., 2., *[d["primary"]["x_match"] for d in branch["detunings"].values()]})
        crossings = [math.acosh(math.sqrt(k)*h/math.sqrt(mu2))/h if k*h*h>mu2 else 0. for k in K2]
        times = sorted(set(times + crossings))
        for x in times:
            for n, f in enumerate(transfer(x, h, mu2)):
                matrices.append([key, x, n, int(K2[n]), int(DEGENERACY[n]), *f.ravel(),
                                 math.cosh(h*x)**-3])
        result["branches"][key] = dict(h=h, mu_squared_over_Htop_squared=mu2,
                                        stationary_theta=top,
                                        initial_strictly_tachyonic_max_n=max(int(n) for n in SHELLS if K2[n]*h*h<mu2),
                                        x_first_nonpositive_frequency_squared_by_n=crossings,
                                        ten_percent_force_radius_minus_plus=radii,
                                        validation=validation(times, h, mu2),
                                        detunings={t:case_result(h,mu2,d,radii)
                                                   for t,d in branch["detunings"].items()})
    destination.mkdir(parents=True, exist_ok=True)
    (destination/"HARMONICS.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    with (destination/"TRANSFER.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["branch", "x", "n", "k_squared", "degeneracy", "q_from_q", "q_from_v",
                         "v_from_q", "v_from_v", "exact_Wronskian"])
        writer.writerows(matrices)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.prior, args.output_dir)
