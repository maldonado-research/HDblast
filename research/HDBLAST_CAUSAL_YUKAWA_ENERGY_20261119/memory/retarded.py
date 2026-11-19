"""Check the twice-subtracted retarded Dirac pair kernel on a native scalar mode."""
import argparse
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import cumulative_simpson, quad, simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline, PchipInterpolator

NORMALIZATION = 8 * np.pi**2
SWITCH = (1.35, 2.05)
READOUTS = (1.7, 2.07, 2.1)


def susceptibility(p2, mu2, max_mass2=None):
    """Value and invariant slope vanish at p²=−mu². Continue from Im ω>0."""
    p2 = np.asarray(p2, dtype=complex)
    safe = np.where(p2 == 0, -mu2, p2)
    if max_mass2 is None:
        return (p2 * (np.log(mu2) - np.log(-safe)) + p2 + mu2) / NORMALIZATION
    log_pivot = np.log1p(max_mass2 / mu2)
    pivot_slope = log_pivot - max_mass2 / (max_mass2 + mu2)
    return (p2 * np.log1p(-max_mass2 / safe) + mu2 * log_pivot
            - (p2 + mu2) * pivot_slope) / NORMALIZATION


def dispersion_quadrature(p2, mu2, top):
    def kernel(v):
        return ((p2 + mu2)**2 * v
                / (NORMALIZATION * (v + mu2)**2 * (v - p2)))
    return (quad(lambda v: kernel(v).real, 0, top, epsabs=2e-12)[0]
            + 1j * quad(lambda v: kernel(v).imag, 0, top, epsabs=2e-12)[0])


class History:
    """Free complex mode on the released fine retained quiet wall in proper time."""

    def __init__(self, native, q_index=None):
        fork = json.loads((native / "FORK.json").read_text())["fine"]
        if q_index is None:
            q_index = len(fork["q"]) // 2
        self.q = fork["q"][q_index]
        with zipfile.ZipFile(native / "WALLS.zip") as archive:
            with archive.open("fine/retain/quiet/wall.csv") as stream:
                wall = np.genfromtxt(stream, names=True, delimiter=",")
        times = wall["H0tau"]
        self.loga = CubicHermiteSpline(times, wall["A_b"] - fork["reference_A"],
                                       wall["H_over_H0"])
        mass = CubicHermiteSpline(times, 100 * (wall["phi_b"] - .5),
                                  100 * wall["v_over_H0"])
        start = float(times[0])
        proper = np.linspace(start, READOUTS[-1], 16385)
        conformal = cumulative_simpson(np.exp(-self.loga(proper)), x=proper, initial=0.)
        self.eta = PchipInterpolator(proper, conformal)
        self.proper = PchipInterpolator(conformal, proper)
        self.lo, self.hi = self.eta(SWITCH)
        inertia = fork["a"]**3 * np.hypot(self.q / fork["a"], mass(start))
        chi0 = 1 / np.sqrt(2 * inertia)

        def rhs(s, state):
            a = np.exp(self.loga(s))
            freq2 = (self.q / a)**2 + mass(s)**2
            return [state[1] / a**3, -a**3 * freq2 * state[0]]

        sol = solve_ivp(rhs, (start, READOUTS[-1]), [chi0 + 0j, -1j * inertia * chi0],
                        dense_output=True, method="DOP853", rtol=2e-11, atol=1e-13,
                        max_step=.0025)
        if not sol.success:
            raise RuntimeError(sol.message)
        self.free = sol.sol

    def source(self, eta):
        t = np.atleast_1d(eta)
        value = np.zeros(t.shape, dtype=complex)
        inside = (self.lo < t) & (t < self.hi)
        proper = self.proper(t[inside])
        if len(proper):
            angle = np.pi * (proper - SWITCH[0]) / (SWITCH[1] - SWITCH[0])
            value[inside] = np.sin(angle)**4 * np.exp(self.loga(proper)) * self.free(proper)[0]
        return value if np.ndim(eta) else value[0]

    def source_derivative(self, t):
        derivative = np.zeros(t.shape, dtype=complex)
        inside = (self.lo < t) & (t < self.hi)
        proper = self.proper(t[inside])
        rate = np.pi / (SWITCH[1] - SWITCH[0])
        angle = rate * (proper - SWITCH[0])
        switch = np.sin(angle)**4
        dswitch = 4 * rate * np.sin(angle)**3 * np.cos(angle)
        a = np.exp(self.loga(proper))
        chi, momentum = self.free(proper)
        derivative[inside] = (a**2 * chi * (dswitch + switch * self.loga(proper, 1))
                              + switch * momentum / a)
        return derivative


def bromwich(history, times, samples, mu2, cutoff=None):
    # Our forward convention is +iωt. The imaginary shift suppresses
    # contamination from the periodically repeated source on the preceding period.
    duration = float(history.eta(READOUTS[-1]))
    dt = duration / samples
    length = 8 * samples
    gamma = 4 / duration
    t = np.arange(length) * dt
    # Only the finite support enters the forward transform.
    damped = history.source(t) * np.exp(-gamma * t)
    frequencies = 2 * np.pi * np.fft.fftfreq(length, dt)
    p2 = (frequencies + 1j * gamma)**2 - history.q**2
    top = None if cutoff is None else cutoff**2 - history.q**2
    response = np.fft.fft(np.fft.ifft(damped) * susceptibility(p2, mu2, top))
    # Undo the Laplace damping only on the requested period of interest.
    stop = samples + 2
    sol = CubicSpline(t[:stop], response[:stop] * np.exp(gamma * t[:stop]))
    return sol(times)


def post_switch_integral(history, time):
    if time <= history.hi:
        raise ValueError("This form of the retarded kernel requires a positive time gap.")

    def integrand(t):
        gap = time - t
        phase = history.q * gap
        kernel = -(np.cos(phase) + phase * np.sin(phase)) / (2 * np.pi**2 * gap**3)
        return kernel * history.source(t)

    return (quad(lambda t: integrand(t).real, history.lo, history.hi, epsabs=2e-10)[0]
            + 1j * quad(lambda t: integrand(t).imag, history.lo, history.hi, epsabs=2e-10)[0])


def direct_pair_energy(history, times, cutoff):
    """Conformal work from the positive spectrum for F = Re(f a χ)."""
    # A uniform Simpson rule in time lets each frequency integral reuse its weights.
    dt = times[1] - times[0]
    quadrature = np.ones(len(times))
    quadrature[1:-1:2] = 4
    quadrature[2:-1:2] = 2
    source = quadrature * (dt / 3) * history.source(times).real
    nodes, weights = np.polynomial.legendre.leggauss(16)
    edges = np.linspace(history.q, cutoff, int(np.ceil((cutoff-history.q)/12))+1)
    total = 0.
    for lo, hi in zip(edges[:-1], edges[1:]):
        omega = (hi+lo)/2 + (hi-lo)/2 * nodes
        amplitude = np.exp(1j * omega[:, None] * times) @ source
        density = omega * (omega**2 - history.q**2) / NORMALIZATION
        total += ((hi-lo)/2) * np.dot(weights, density * abs(amplitude)**2)
    return float(total)


def work_sum_rule(history, mu2):
    times = np.linspace(history.lo, history.hi, 4097)
    derivative = history.source_derivative(times).real

    def response_work(pivot, cutoff):
        response = bromwich(history, times, 4096, pivot, cutoff)
        return float(-simpson(derivative * response.real, x=times))

    rows = {}
    for cutoff in (40,80,160,400,800):
        rows[str(cutoff)] = {
            "positive_pair_spectrum": direct_pair_energy(history, times, cutoff),
            "retarded_work": response_work(mu2, cutoff),
        }
    rows["infinite"] = {
        "retarded_work": response_work(mu2, None),
        "retarded_work_at_mu_4": response_work(4*mu2, None),
    }
    return rows


def pair(z):
    return [float(z.real), float(z.imag)]


def run(native):
    history = History(native)
    times = history.eta(READOUTS)
    mu = 2.
    checks = {}
    for p2 in ((2+1j)**2-history.q**2, (9+.8j)**2-history.q**2):
        for top in (1600., np.inf):
            calculated = susceptibility(p2, mu**2, None if np.isinf(top) else top)
            other = dispersion_quadrature(p2, mu**2, top)
            checks["dispersion_max_absolute_error"] = max(
                checks.get("dispersion_max_absolute_error", 0.), float(abs(other - calculated)))
    reference = bromwich(history, times, 4096, mu**2)
    coarse = bromwich(history, times, 2048, mu**2)
    checks["grid_refinement_max_absolute_error"] = float(np.max(abs(reference - coarse)))
    alternative = bromwich(history, times, 4096, (2*mu)**2)
    tail = [post_switch_integral(history, t) for t in times[1:]]
    checks["post_switch_quadrature_max_absolute_error"] = float(max(abs(reference[1:] - tail)))
    checks["post_switch_scale_change_max_absolute"] = float(max(abs(alternative[1:] - reference[1:])))
    mid = times[0]
    h = 2e-4
    f = history.source(mid)
    side = history.source(mid + h * np.array([-2,-1,1,2]))
    d2 = (-side[0]+16*side[1]-30*f+16*side[2]-side[3]) / (12*h**2)
    predicted_shift = (np.log(4.) * (-d2 - history.q**2 * f) + 3*mu**2 * f) / NORMALIZATION
    checks["interior_local_scale_shift_error"] = float(abs(alternative[0]-reference[0]-predicted_shift))
    acausal = bromwich(history, history.eta([1.31, 1.32]), 4096, mu**2)
    checks["before_switch_max_absolute"] = float(max(abs(acausal)))
    # Positive retarded cut at ω>q and negative-frequency reflection symmetry.
    omega = 9.
    epsilon = 1e-9
    expected_imag = np.pi * (omega**2-history.q**2) / NORMALIZATION
    positive = susceptibility((omega+1j*epsilon)**2-history.q**2, mu**2)
    negative = susceptibility((-omega+1j*epsilon)**2-history.q**2, mu**2)
    checks["imaginary_cut_absolute_error"] = float(abs(positive.imag - expected_imag))
    checks["negative_frequency_conjugation_error"] = float(abs(negative-np.conj(positive)))
    checks["light_cone_absolute_error"] = float(abs(susceptibility(0., mu**2) - mu**2/NORMALIZATION))
    work = work_sum_rule(history, mu**2)
    checks["finite_work_identity_max_absolute_error"] = max(
        abs(row["positive_pair_spectrum"] - row["retarded_work"])
        for key,row in work.items() if key != "infinite")
    checks["infinite_work_vs_high_frequency_reference"] = abs(
        work["infinite"]["retarded_work"]-work["800"]["positive_pair_spectrum"])
    checks["infinite_work_scale_change_max_absolute"] = abs(
        work["infinite"]["retarded_work_at_mu_4"]-work["infinite"]["retarded_work"])
    if (checks["dispersion_max_absolute_error"] > 2e-10
            or max(checks.values()) > 1e-5):
        raise ArithmeticError(checks)
    finite = {}
    for cutoff in (40,80,160,320,640):
        result = bromwich(history, times, 4096, mu**2, cutoff)
        finite[str(cutoff)] = {
            "response": [pair(z) for z in result],
            "absolute_error_from_infinite_cutoff": [float(x) for x in abs(result-reference)],
        }
    return {
        "wall": "fine/retain/quiet", "native_q": history.q, "switch_proper": SWITCH,
        "mu": mu, "proper_readout": READOUTS, "conformal_readout_since_fork": list(times),
        "infinite_cutoff_response": [pair(z) for z in reference],
        "infinite_cutoff_response_at_mu_4": [pair(z) for z in alternative],
        "post_switch_direct_time_integral": [pair(z) for z in tail],
        "finite_comoving_energy_cutoffs": finite,
        "conformal_dissipated_work_real_mode": work,
        "checks": checks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=Path(__file__).resolve().parents[2]
                        / "HDBLAST_COVARIANCE_FEEDBACK_20261119" / "native")
    print(json.dumps(run(parser.parse_args().native), indent=2, allow_nan=False))
