"""First-collision upper bounds for isotropic massless particle–antiparticle scattering."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre


MEMORY = Path("research/HDBLAST_QUANTUM_MEMORY_20261119/quantum")
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
ARMS = ("quiet", "plus", "minus", "far")
READOUTS = (2.10, 2.16)


def load_spectrum(repo):
    path = repo / MEMORY / "spectrum.py"
    spec = importlib.util.spec_from_file_location("archived_ordered_spectrum", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rule(count):
    x, w = roots_legendre(count)
    return (x + 1) / 2, w / 2


def moments(fork, profile, band, order, squared_threshold):
    """Integrate the ordered pair yield, analytically averaging inverse daughter energy."""
    q = np.asarray(fork["q"])
    wq = np.asarray(fork["weight"])
    eta, dt, h, k = profile
    u, du = rule(order)
    X, C, P = (np.asarray(fork["states"]["retain"][v]) for v in ("X", "C", "P"))
    summed = np.zeros(3)
    for j, momentum in enumerate(q):
        lower, upper = (momentum, 60.) if band == "low" else (60., 120.)
        width = upper - lower
        if squared_threshold and band == "low":
            omega, dw = momentum + width * u**2, 2 * width * u * du
        else:
            omega, dw = lower + width * u, width * du
        amp = np.exp(1j * omega[:, None] * eta) @ np.stack((dt*h[j], dt*k[j]), axis=1)
        v, z = amp.T
        # Completing a square implements v† W v with W_12=C+i/2.
        ordered = X[j] * abs(v + (C[j]+0.5j)*z/X[j])**2
        ordered += (P[j]-(C[j]**2+0.25)/X[j]) * abs(z)**2
        if np.min(ordered) < -2e-15:
            raise ValueError("Negative ordered spectral power")
        pair = wq[j] * dw * ((omega-momentum)*(omega+momentum)/(8*np.pi**2)) * ordered
        inverse = np.log1p(2*momentum/(omega-momentum)) / momentum
        summed += [np.sum(pair), np.sum(pair * omega / 2), np.sum(pair * inverse)]
    return summed


def relative(a, b):
    return float(max(abs(a-b)/abs(b)))


def readout(m, column, a, h):
    n, e, inv = m
    # Number-weighted and energy-weighted averages over the one-particle marginal.
    tau_n = float(8*np.pi * column * inv**2 / n)
    tau_e = float(8*np.pi * column * inv*n / e)
    return dict(
        mean_tau_number_per_y2=tau_n, mean_tau_energy_per_y2=tau_e,
        necessary_y2_for_unit_mean_tau=1./tau_n,
        frac_number_at_y_0_1_upper=float(-np.expm1(-0.01*tau_n)),
        frac_energy_at_y_0_1_upper=float(-np.expm1(-0.01*tau_e)),
        mean_rate_over_H_number_per_y2=float(8*np.pi*inv**2/(n*a*h)),
        mean_rate_over_H_energy_per_y2=float(8*np.pi*inv*n/(e*a*h)),
    )


def evolving_targets(m, column):
    """Upper fraction of original tags ever hit, with arbitrary elastic redistribution."""
    y = 0.1
    g = 2.
    n, e, initial_inverse = y**2 * m
    fermi_momentum = (6*np.pi**2*n/g)**(1./3.)
    pauli_inverse = 1.5*n/fermi_momentum
    coefficient = pauli_inverse/initial_inverse
    exposure = 8*np.pi * column*initial_inverse**2/n

    def primitive(fraction):
        root = np.cbrt(fraction)
        return 3*(root - np.arctan(np.sqrt(coefficient)*root)/np.sqrt(coefficient))/coefficient

    upper = brentq(lambda f: primitive(f)-exposure, 0., 1., xtol=5e-18)
    solution = solve_ivp(lambda _, f: 1 + coefficient*np.maximum(f,0.)**(2./3.),
                         (0.,exposure),[0.],method="DOP853",atol=1e-17,rtol=1e-12,
                         max_step=exposure/100.)
    if not solution.success or abs(solution.y[0,-1]-upper)>2e-13:
        raise ValueError("Implicit bound and comparison differential equation disagree")
    # Quadrature of the positive inverse derivative also checks the subtraction in primitive.
    check = quad(lambda t: 3*t*t/(1+coefficient*upper**(2./3.)*t*t),0.,1.)[0] * upper
    if abs(check-exposure)>2e-17:
        raise ValueError("Integral of reciprocal rate")
    min_column_ten_percent = n * primitive(0.1) / (8*np.pi*initial_inverse**2)
    agreement = quad(lambda f: 1/(1+coefficient*f**(2/3)),0.,0.1,epsabs=5e-14)[0]
    if abs(agreement-primitive(0.1))>5e-13:
        raise ValueError("Necessary ten-percent column integral")
    energy_depth = upper*n**2/(e*initial_inverse)
    any_target_depth = 8*np.pi*column*pauli_inverse*initial_inverse/n
    return dict(production_y=y, spin_states_per_charge=int(g), density_per_charge=n,
                initial_inverse_moment=initial_inverse, pauli_inverse_moment=pauli_inverse,
                comparison_coefficient=coefficient, initial_exposure=exposure,
                minimum_conformal_column_for_ten_percent=float(min_column_ten_percent),
                minimum_column_over_observed_for_ten_percent=float(min_column_ten_percent/column),
                remaining_column_necessary_for_ten_percent=float(max(0.,min_column_ten_percent-column)),
                first_collision_fraction_number_upper=float(upper),
                first_collision_fraction_initial_energy_upper=float(-np.expm1(-energy_depth)),
                all_density_pauli_number_upper=float(-np.expm1(-any_target_depth)),
                independent_ode_discrepancy=float(abs(solution.y[0,-1]-upper)))


def experiment(repo):
    spectrum = load_spectrum(repo)
    previous = json.loads((repo / MEMORY / "SPECTRUM.json").read_text())
    fork = json.loads((repo / NATIVE / "FORK.json").read_text())["fine"]
    quiet = spectrum.wall(repo / NATIVE, "quiet", fork)[1]
    qmax = max(fork["q"])
    out = {}
    for arm in ARMS:
        loga, phi = spectrum.wall(repo / NATIVE, arm, fork)
        ends, solution, _ = spectrum.fundamental(fork, loga, phi, quiet)
        start = ends[1]
        fine = spectrum.profile(fork, ends, solution, phi, quiet, 5)
        coarse = spectrum.profile(fork, ends, solution, phi, quiet, 3)
        feeds = {}
        controls = {}
        for band in ("low", "high"):
            feed = moments(fork, fine, band, 144, True)
            refined = moments(fork, coarse, band, 96, True)
            plain = moments(fork, fine, band, 192, False)
            controls[band] = dict(relative_change_under_refinement=relative(refined,feed),
                                  ordinary_Gauss_relative_difference=relative(plain,feed))
            if relative(refined, feed) > 5e-8 or relative(plain, feed) > 8e-6:
                raise ValueError(f"Spectral convergence for {arm}/{band}: {controls[band]}")
            feeds[band] = feed
        feeds["all"] = feeds["low"] + feeds["high"]
        native_ref = previous["arms"][arm]
        for band, cut in (("low", "60"), ("all", "120")):
            energy = 2*feeds[band][1]/native_ref["scale_factor_at_readout"]
            recorded = native_ref["positive_frequency_cutoff_checks"][cut]["archived"]
            if abs(energy - recorded) > 6e-7:
                raise ValueError("Independent archived energy check")
        for n, e, inv in feeds.values():
            if not (n > 0 and inv*e >= n*n):
                raise ValueError("Cauchy–Schwarz moment bound")
        n0, _, inv0 = feeds["high"]
        if not n0/(60.+qmax/2) <= inv0 <= n0/(30.-qmax/2):
            raise ValueError("High-band support")
        moments_out = {b:dict(zip(("N","M1","Mminus1"),map(float,m))) for b,m in feeds.items()}
        reads = {}
        for s in READOUTS:
            inside = loga.x[(loga.x>start) & (loga.x<s)]
            spans = np.r_[start, inside, s]
            column = sum(quad(lambda t: np.exp(-loga(t)), l, r)[0] for l,r in zip(spans[:-1],spans[1:]))
            u, w = rule(8)
            lengths = np.diff(spans)
            points = spans[:-1,None] + lengths[:,None]*u
            check = float(np.sum(lengths[:,None]*w*np.exp(-loga(points))))
            if abs(column-check) > 5e-14:
                raise ValueError("Conformal time integral")
            a,h = float(np.exp(loga(s))),float(loga(s,1))
            # The high subset sees all targets in high_on_all, as opposed to internal high/high.
            diagnostics = {band:readout(m,column,a,h) for band,m in feeds.items()}
            for tag in ("low","high"):
                n,e,inv=feeds[tag]
                target=feeds["all"][2]
                diagnostics[tag+"_on_all"]=dict(mean_tau_number_per_y2=float(8*np.pi*column*inv*target/n),
                                                 mean_tau_energy_per_y2=float(8*np.pi*column*n*target/e))
            frac=feeds["high"][2]/feeds["all"][2]
            shares=dict(low_low=float((1-frac)**2),low_high_and_high_low=float(2*frac*(1-frac)),
                        high_high=float(frac**2))
            if abs(sum(shares.values())-1)>2e-15:
                raise ValueError("Channel accounting")
            reads[str(s)] = dict(conformal_path=float(column),scale_factor=a,expansion_rate_h=h,
                                 rates=diagnostics,shares_of_number_weighted_envelope=shares,
                                 evolving_targets=evolving_targets(feeds["all"],column))
        out[arm]=dict(gate_closes=float(start),moments_per_species=moments_out,convergence=controls,
                      readouts=reads)
    return dict(grid="fine",target="one massless Dirac antiparticle species, spin summed",
                cross_section_assumption="sigma_tot = sigma_tr <= 16 pi / s_hat; isotropic elastic single-wave envelope",
                histories=out)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(experiment(args.repo),indent=2,allow_nan=False)+"\n")
