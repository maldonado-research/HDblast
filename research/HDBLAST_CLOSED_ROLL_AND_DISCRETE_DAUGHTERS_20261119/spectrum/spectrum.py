#!/usr/bin/env python3
"""Finite, smooth homogeneous pulses producing a scalar on an externally held round S3.

Units inside the mode solver are m_sigma=1. Radius a is quoted in H_top^-1,
so its dimensionless solver value is a*mu, mu=m_sigma/H_top. The static sphere
and prescribed pulse are a spectral benchmark, not an Einstein-equation solution.
"""
import argparse
import csv
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import roots_legendre


def packet(t, cycles, power=4):
    """f=sin^power(pi*t/T)cos(t), T=2*pi*N; used only within [0,T]."""
    h = .5 / cycles
    s, c = np.sin(h*t), np.cos(h*t)
    w = s**power
    wp = h*power*s**(power-1)*c
    wpp = h*h*power*((power-1)*s**(power-2)*c*c-s**power)
    carrier, carrier_prime = np.cos(t), -np.sin(t)
    return w*carrier, wp*carrier+w*carrier_prime, wpp*carrier+2*wp*carrier_prime-w*carrier


def effective_time(cycles, power=4):
    return 2*np.pi*cycles*comb(2*power, power)/4**power


def born_beta(k, cycles, epsilon=.001, power=4):
    """Exact compact-pulse Fourier transform, including the negative carrier."""
    k = np.asarray(k, dtype=float)
    T = 2*np.pi*cycles
    transform = np.zeros(k.shape, complex)
    for j in range(power+1):
        offset = (j-power/2)/cycles
        coefficient = (-1)**(j-power//2)*comb(power, j)/2**power
        for carrier in (-1, 1):
            mismatch = 2*k-carrier-offset
            transform += (.5*T*coefficient*np.exp(-.5j*T*mismatch)
                          *np.sinc(.5*T*mismatch/np.pi))
    return 2j*epsilon*k*transform


def exact_modes(k, cycles, epsilon=.001, power=4, rtol=5e-11, atol=3e-14):
    """Bogoliubov equations relative to static in/out frequency k>0."""
    k = np.asarray(k, dtype=float)
    if np.any(k <= 0):
        raise ValueError("The massless zero mode has no oscillator vacuum")
    initial = np.r_[np.ones(k.size, complex), np.zeros(k.size, complex)]
    def rhs(t, y):
        _, fp, fpp = packet(t, cycles, power)
        V = epsilon*fpp + (epsilon*fp)**2
        alpha, beta = np.split(y, 2)
        wave = np.exp(-2j*k*t)
        source = V/(2j*k)
        return np.r_[-source*(alpha+wave.conj()*beta), source*(wave*alpha+beta)]
    sol = solve_ivp(rhs, (0., 2*np.pi*cycles), initial, method="DOP853", t_eval=[2*np.pi*cycles],
                    rtol=rtol, atol=atol, max_step=.7)
    if not sol.success:
        raise RuntimeError(sol.message)
    return np.split(sol.y[:, -1], 2)


def direct_beta(k, cycles, epsilon=.001, power=4):
    """Independent second-order oscillator, normalized incoming mode exp(-ikt)."""
    k = np.asarray(k, dtype=float)
    def rhs(t, y):
        q, v = np.split(y, 2)
        _, fp, fpp = packet(t, cycles, power)
        V = epsilon*fpp+(epsilon*fp)**2
        return np.r_[v, (V-k*k)*q]
    T = 2*np.pi*cycles
    sol = solve_ivp(rhs, (0, T), np.r_[np.ones(k.size, complex), -1j*k],
                    method="DOP853", t_eval=[T], rtol=2e-12, atol=2e-14, max_step=.2)
    if not sol.success:
        raise RuntimeError(sol.message)
    q, v = np.split(sol.y[:, -1], 2)
    return .5*np.exp(-1j*k*T)*(q+v/(1j*k))


def continuum_grid(cycles, kmax=4., order=8):
    # Each shifted sinc oscillation is resolved independently.
    edges = np.arange(round(2*cycles*kmax)+1, dtype=float)/(2*cycles)
    u, w = roots_legendre(order)
    half = np.diff(edges)[:, None]/2
    return (edges[:-1, None]+half*(1+u)).ravel(), (half*w).ravel()


def densities(k, occupation, weights):
    number = float(weights@occupation)
    energy = float((weights*k)@occupation)
    return dict(number_density_over_m3=number, energy_density_over_m4=energy,
                mean_energy_over_m=energy/number if number else None,
                maximum_occupation=float(np.max(occupation)),
                pair_density_over_m3=number/2)


def continuum(cycles, epsilon=.001, power=4, order=8, kmax=4., solve_exact=True):
    k, dk = continuum_grid(cycles, kmax, order)
    weights = dk*k*k/(2*np.pi**2)
    born = densities(k, abs(born_beta(k, cycles, epsilon, power))**2, weights)
    result = dict(born=born, cutoff_over_m=kmax, quadrature_order=order)
    if solve_exact:
        alpha, beta = exact_modes(k, cycles, epsilon, power)
        result["exact"] = densities(k, abs(beta)**2, weights)
        result["maximum_Wronskian_residual"] = float(np.max(abs(abs(alpha)**2-abs(beta)**2-1)))
    # This first-order number identity is independent of the frequency quadrature.
    T = 2*np.pi*cycles
    time_int = quad(lambda t: (epsilon*packet(t, cycles, power)[2])**2,
                    0., T, epsabs=1e-17, epsrel=2e-12, limit=500)[0]/(16*np.pi)
    result["born_number_Parseval_full_line_over_m3"] = time_int
    result["born_number_Parseval_fractional_difference"] = born["number_density_over_m3"]/time_int-1
    if abs(result["born_number_Parseval_fractional_difference"]) > 1e-7:
        raise ValueError("Continuum particle-number Parseval identity failed")
    result["flat_width_energy_over_m4"] = epsilon**2*effective_time(cycles, power)/(64*np.pi)
    return result


def spectrum_grid(mu, radius, kmax=4.):
    r = mu*radius
    nmax = int(np.floor(np.sqrt((r*kmax)**2+1)-1))
    n = np.arange(1, nmax+1)
    return n, np.sqrt(n*(n+2))/r, (n+1.)**2/(2*np.pi**2*r**3)


def born_energy_tail_bound(cycles, radius_m, cutoff=4., epsilon=.001):
    """Absolute ultraviolet bounds for p=4, from its rational Fourier transform.

    The bounds apply to Born energy only; exact-mode cutoff checks are separate.
    j=n+1 and j0 is the first omitted spatial level. Both bounds are in m^4.
    """
    b = 1/cycles
    nu = 2*cutoff
    if nu <= 1+2*b:
        raise ValueError("Tail estimate requires frequencies above all pulse harmonics")
    B = 3*b**4/((1-1/nu)**5*(1-(b/(nu-1))**2)*(1-(2*b/(nu-1))**2))
    j0 = int(np.floor(np.sqrt(1+(radius_m*cutoff)**2)))+1
    continuum_bound = epsilon**2*B**2/(128*np.pi**2*nu**4)
    sphere_bound = (epsilon**2*B**2*radius_m**4/(512*np.pi**2)
                    *(1-1/j0**2)**(-3.5)*(j0**-5+.25*j0**-4))
    return dict(continuum_energy_over_m4=continuum_bound,
                discrete_energy_over_m4=sphere_bound, first_omitted_j=j0)


def calculate_case(mu, radius, cycles, epsilon=.001, power=4, kmax=4., solve_exact=True):
    n, k, weights = spectrum_grid(mu, radius, kmax)
    born_occ = abs(born_beta(k, cycles, epsilon, power))**2
    result = dict(radius_Htop_inverse=radius, radius_m_inverse=mu*radius,
                  volume_Htop_inverse_cubed=2*np.pi**2*radius**3,
                  volume_m_inverse_cubed=2*np.pi**2*(mu*radius)**3,
                  cycles=cycles, epsilon=epsilon, window_sine_power=power,
                  duration_Htop_inverse=2*np.pi*cycles/mu, duration_m_inverse=2*np.pi*cycles,
                  effective_duration_Htop_inverse=effective_time(cycles, power)/mu,
                  cutoff_over_m=kmax, n_max=int(n[-1]),
                  born=densities(k, born_occ, weights))
    nearest = int(np.argmin(abs(2*k-1)))
    result["nearest_level"] = dict(n=int(n[nearest]), degeneracy=int((n[nearest]+1)**2),
                                 omega_over_m=float(k[nearest]),
                                 fractional_pair_detuning=float(2*k[nearest]-1),
                                 detuning_times_duration=float((2*k[nearest]-1)*2*np.pi*cycles))
    exact_occ = np.full(k.size, np.nan)
    if solve_exact:
        alpha, beta = exact_modes(k, cycles, epsilon, power)
        exact_occ = abs(beta)**2
        result["exact"] = densities(k, exact_occ, weights)
        result["maximum_Wronskian_residual"] = float(np.max(abs(abs(alpha)**2-abs(beta)**2-1)))
        result["exact_to_Born_energy"] = result["exact"]["energy_density_over_m4"]/result["born"]["energy_density_over_m4"]
    rows = [dict(n=int(ni), degeneracy=int((ni+1)**2), omega_over_m=float(ki),
                 born_occupation=float(bi), exact_occupation=float(ei))
            for ni, ki, bi, ei in zip(n, k, born_occ, exact_occ)]
    return result, rows


def load_inputs(gravity):
    gravity = Path(gravity).resolve()
    ledger_path = gravity/"previous/prior/archive/BLAST_REHEATING_RESULTS.json"
    background_path = gravity.parent/"modulus/background/BACKGROUND.json"
    widths_path = gravity.parent/"modulus/widths/WIDTHS.json"
    ledger = json.loads(ledger_path.read_text())["B"]
    background = json.loads(background_path.read_text())["branches"]
    widths = json.loads(widths_path.read_text())["branches"]
    result = {}
    for label in ("phi_m=0.5", "phi_m=-0.5"):
        b = ledger[label]
        mu = b["m_over_H_top"]
        alpha = b["alpha"]
        if mu != widths[label]["m_over_Htop"] or alpha != widths[label]["alpha"]:
            raise ValueError("Archival mass/trace-coupling mismatch")
        closed = background[label]["closed"]["off"]
        result[label] = dict(phi_minimum=b["phi_m"], mu=mu, alpha=alpha,
                             canonical_amplitude_Mpl=.001/alpha,
                             published_minimum_crossing_radius=closed["first_minimum_crossing"]["a"],
                             published_first_maximum_radius=closed["turnaround"]["a"])
    return result


def experiment(gravity):
    branches = load_inputs(gravity)
    continua = {str(N): continuum(N) for N in (4, 8, 16)}
    cases, rows = [], []
    for label, b in branches.items():
        mu = b["mu"]
        # Pick the resonant level near a=1 independently, then construct exact
        # resonance and a midpoint *in frequency*, not in radius.
        level = max(1, round(np.sqrt(1+(mu/2)**2)-1))
        e0, e1 = np.sqrt(level*(level+2)), np.sqrt((level+1)*(level+3))
        radii = dict(unit=1., representative=1.1,
                     published_minimum_crossing=b["published_minimum_crossing_radius"],
                     published_first_maximum=b["published_first_maximum_radius"],
                     matched_resonance=2*e0/mu, matched_midpoint=(e0+e1)/mu)
        b["matched_control_level"] = level
        for role, radius in radii.items():
            for N in (4, 8, 16):
                case, mode_rows = calculate_case(mu, radius, N)
                case.update(branch=label, radius_role=role)
                ref = continua[str(N)]
                for kind in ("born", "exact"):
                    for moment, key in (("number", "number_density_over_m3"), ("energy", "energy_density_over_m4")):
                        case[f"ratio_{moment}_{kind}_to_flat_continuum"] = case[kind][key]/ref[kind][key]
                    case[kind]["number_density_over_Htop3"] = case[kind]["number_density_over_m3"]*mu**3
                    case[kind]["energy_density_over_Htop4"] = case[kind]["energy_density_over_m4"]*mu**4
                case["ratio_exact_energy_to_flat_width"] = case["exact"]["energy_density_over_m4"]/ref["flat_width_energy_over_m4"]
                cases.append(case)
                for row in mode_rows:
                    rows.append(dict(branch=label, radius_role=role, radius_Htop_inverse=radius, cycles=N, **row))
    controls = {}
    # All continuum controls are universal in m_sigma units.
    refined = continuum(16, order=12, solve_exact=False)
    ultraviolet = {str(kmax): continuum(4, order=12, kmax=kmax, solve_exact=False)["born"]
                   for kmax in (4., 8., 16.)}
    controls["continuum_refined_N16"] = refined
    controls["continuum_Born_N4_by_cutoff"] = ultraviolet
    exact_ultraviolet = continuum(4, order=8, kmax=8.)
    controls["continuum_exact_N4_cutoff8"] = exact_ultraviolet
    controls["continuum_exact_N4_cutoff4_to8_energy_fractional_change"] = (
        exact_ultraviolet["exact"]["energy_density_over_m4"]/continua["4"]["exact"]["energy_density_over_m4"]-1)
    controls["maximum_continuum_quadrature_fractional_change"] = max(
        abs(refined["born"][key]/continua["16"]["born"][key]-1)
        for key in ("number_density_over_m3", "energy_density_over_m4"))
    controls["discrete_amplitude_and_sign"] = {}
    controls["discrete_Born_by_cutoff"] = {}
    controls["large_radius_Born"] = {}
    controls["analytic_resonance_occupation"] = {}
    controls["Born_energy_above_cutoff4_bounds"] = {}
    controls["discrete_exact_N4_cutoff8"] = {}
    direct_errors, null_max, sign_errors = [], [], []
    for label, b in branches.items():
        mu = b["mu"]
        n, k, weights = spectrum_grid(mu, 1.1)
        amplitude = {}
        for eps in (0., .0005, .001, .002, -.001):
            a, beta = exact_modes(k, 8, eps)
            amplitude[str(eps)] = densities(k, abs(beta)**2, weights)
            if eps == 0.:
                null_max.append(float(np.max(abs(beta))))
            if eps == .001:
                positive = a, beta
            if eps == -.001:
                sign_errors += [float(np.max(abs(a-positive[0]))), float(np.max(abs(beta+positive[1])))]
        controls["discrete_amplitude_and_sign"][label] = amplitude
        # Direct oscillator comparison at several low levels and the resonant neighborhood.
        take = np.unique(np.r_[0, 1, np.argmin(abs(k-.5)), min(len(k)-1, 12)])
        _, beta = exact_modes(k[take], 8)
        direct_errors.append(float(np.max(abs(beta-direct_beta(k[take], 8)))))
        controls["discrete_Born_by_cutoff"][label] = {
            str(cut): calculate_case(mu, 1.1, 4, kmax=cut, solve_exact=False)[0]["born"]
            for cut in (4., 8., 16.)}
        controls["discrete_exact_N4_cutoff8"][label] = calculate_case(mu, 1.1, 4, kmax=8.)[0]
        controls["large_radius_Born"][label] = {
            str(rad): {key: calculate_case(mu, rad, 16, solve_exact=False)[0]["born"][key]/continua["16"]["born"][key]
                       for key in ("number_density_over_m3", "energy_density_over_m4")}
            for rad in (10., 20., 40.)}
        controls["Born_energy_above_cutoff4_bounds"][label] = {
            str(N): born_energy_tail_bound(N, 1.1*mu) for N in (4, 8, 16)}
    for N in (4, 8, 16):
        analytic = (3*np.pi*.001*N/8)**2
        transform_value = float(abs(born_beta(np.array([.5]), N)[0])**2)
        controls["analytic_resonance_occupation"][str(N)] = dict(
            expected=analytic, analytic_transform=transform_value,
            fractional_difference=transform_value/analytic-1)
    controls["maximum_direct_complex_beta_error"] = max(direct_errors)
    controls["zero_drive_maximum_abs_beta"] = max(null_max)
    controls["reciprocal_metric_maximum_complex_coefficient_error"] = max(sign_errors)
    # An independent time-domain check of the analytic compact Fourier transform.
    transform_errors = []
    for N in (4, 8, 16):
        T = 2*np.pi*N
        for k in (.17, .5, .63, 1.2):
            integral = sum(quad(lambda t: packet(t, N)[0]*f(2*k*t), 0, T,
                                epsabs=2e-11, epsrel=2e-11, limit=500)[0]*c
                           for f, c in ((np.cos, 1), (np.sin, -1j)))
            transform_errors.append(abs(born_beta(np.array([k]), N)[0]-2j*.001*k*integral))
    controls["maximum_analytic_Born_beta_vs_time_quadrature_abs_error"] = float(max(transform_errors))
    maximum_wronskian = max([r["maximum_Wronskian_residual"] for r in cases]
                            +[r["maximum_Wronskian_residual"] for r in continua.values()])
    controls["maximum_Wronskian_residual"] = maximum_wronskian
    if maximum_wronskian > 1e-8 or max(direct_errors) > 2e-8 or max(null_max) != 0. or max(sign_errors) > 2e-8:
        raise ValueError("Oscillator null, phase, or unitary controls failed")
    if controls["maximum_continuum_quadrature_fractional_change"] > 1e-7:
        raise ValueError("Continuum quadrature did not converge")
    if abs(controls["continuum_exact_N4_cutoff4_to8_energy_fractional_change"]) > 1e-7:
        raise ValueError("Exact continuum ultraviolet cutoff did not converge")
    return dict(model="static externally held round S3; one real massless minimally coupled daughter; A/Am=exp(epsilon sin^4(pi tau/T) cos tau)",
                convention="tau=m_sigma t_E; mu=m_sigma/Htop; radius in Htop^-1; occupations count individual particles; n>=1",
                zero_mode="n=0 has omega=0 and no static oscillator ground state; excluded from particle and energy sums; requires separately specified quantum state",
                branches=branches, continua=continua, cases=cases, controls=controls), rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gravity", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result, rows = experiment(args.gravity)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"SPECTRUM.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    with (args.output/"SPECTRUM.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(dict(cases=len(result["cases"]), mode_rows=len(rows), controls=result["controls"]), indent=2))


if __name__ == "__main__":
    main()
