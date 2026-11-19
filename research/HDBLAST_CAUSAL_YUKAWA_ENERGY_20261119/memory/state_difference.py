"""Continuum retarded-work check of state-dependent quantum radiation energy."""
import argparse
import json
import importlib.util
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

spec = importlib.util.spec_from_file_location("retarded", Path(__file__).with_name("retarded.py"))
retarded = importlib.util.module_from_spec(spec)
spec.loader.exec_module(retarded)


def causal_components(history, times, samples, mu):
    derivative = history.source_derivative(times)
    response = retarded.bromwich(history, times, samples, mu**2)
    # Real and imaginary mode quadratures and their mixed symmetric product.
    integrands = (derivative.real * response.real, derivative.imag * response.imag,
                  derivative.real * response.imag + derivative.imag * response.real)
    return -simpson(np.array(integrands), x=times, axis=1)


def spectral_components(history, times, cutoff):
    # Compute two first-order transition amplitudes at each ω independently of
    # the causal response. Simpson integration of a compact time history.
    simpson_weights = np.ones(len(times))
    simpson_weights[1:-1:2] = 4
    simpson_weights[2:-1:2] = 2
    mode = history.source(times)
    dt = times[1] - times[0]
    sources = np.stack((mode, mode.conj()), axis=1) * (dt / 3 * simpson_weights[:,None])
    nodes, weights = np.polynomial.legendre.leggauss(14)
    edges = np.linspace(history.q, cutoff, int(np.ceil((cutoff-history.q)/16))+1)
    total = np.zeros(4)
    for lo, hi in zip(edges[:-1], edges[1:]):
        omega = (hi+lo)/2 + (hi-lo)/2 * nodes
        a, b = (np.exp(1j * omega[:,None] * times) @ sources).T
        density = omega * (omega**2-history.q**2) / retarded.NORMALIZATION
        cross = a * b.conj()
        moments = np.array((abs(a)**2+abs(b)**2, 2*cross.real, -2*cross.imag, abs(b)**2))
        total += moments @ ((hi-lo)/2 * weights * density)
    return total


def run(native):
    fork = json.loads((native / "FORK.json").read_text())["fine"]
    labels = tuple(fork["states"])
    integrated = {k:np.zeros(4) for k in labels}
    vacuum = 0.
    for i, weight in enumerate(fork["weight"]):
        history = retarded.History(native, q_index=i)
        times = np.linspace(history.lo, history.hi, 2049)
        a_final = float(np.exp(history.loga(retarded.READOUTS[-1])))
        measure = weight / a_final
        work = [causal_components(history, times, 2048, mu) for mu in (2., 4.)]
        spectra = [spectral_components(history, times, top) for top in (160., 400.)]
        vacuum += measure * spectra[1][3]
        for label in labels:
            state = fork["states"][label]
            n, re, im = (state[k][i] for k in ("n", "Re_c", "C"))
            quadratures = 2 * np.array((n+re,n-re,-im))
            quantum = np.array((n,re,im))
            differences = ([quadratures @ w for w in work]
                           + [quantum @ spectrum[:3] for spectrum in spectra])
            integrated[label] += measure * np.array(differences)
    checks = {
        "largest_absolute_change_with_scale": max(abs(x[0]-x[1]) for x in integrated.values()),
        "largest_absolute_difference_causal_vs_direct": max(abs(x[0]-x[3]) for x in integrated.values()),
        "largest_absolute_direct_cutoff_difference": max(abs(x[2]-x[3]) for x in integrated.values())
    }
    if max(checks.values()) > 2e-6:
        raise ArithmeticError(checks)
    columns = ("retarded_mu2","retarded_mu4","direct_to_160","direct_to_400")
    return {
        "wall": "fine/retain/quiet", "a_final": a_final,
        "matched_fork_vacuum_energy_direct_to_400": vacuum,
        "state_excess_energy_at_fixed_wall": {k:dict(zip(columns,x)) for k,x in integrated.items()},
        "checks": checks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, default=Path(__file__).resolve().parents[2]
                        / "HDBLAST_COVARIANCE_FEEDBACK_20261119" / "native")
    print(json.dumps(run(parser.parse_args().native), indent=2, allow_nan=False))
