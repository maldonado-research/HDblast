"""Leading switched Yukawa emission from a supplied wall and a fork covariance."""
from pathlib import Path
import hashlib
import importlib.util

import numpy as np
from scipy.integrate import cumulative_simpson, simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

FORK_DRIVER = "research/HDBLAST_MATCHED_PARENT_SOURCES_20261119/source_check/source_transfer.py"
FORK_DIGEST = "1f590e627c8449e40f4500b4434ddc781c72264a7a47d8ab40a94195ead9040b"
WINDOWS = {"wide": (1.35, 2.05), "narrow": (1.45, 1.95)}
CUT = 2.10
CHANNELS = ("vacuum", "diagonal", "real", "imag")


def published_method(repo):
    path = Path(repo)/FORK_DRIVER
    if hashlib.sha256(path.read_bytes()).hexdigest() != FORK_DIGEST:
        raise ValueError("Unexpected preceding source-matching calculation")
    spec = importlib.util.spec_from_file_location("source_matching", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.load_reference(repo)


def wall_interpolants(data, a_reference):
    s = np.asarray(data["H0tau"], float)
    if np.any(np.diff(s) <= 0) or not s[0] < min(lo for lo, _ in WINDOWS.values()) < CUT <= s[-1]:
        raise ValueError("Wall does not cover fork, both switches, and readout")
    loga = CubicHermiteSpline(s, data["A_b"] - a_reference, data["H_over_H0"])
    mass = CubicHermiteSpline(s, 100*(data["phi_b"]-.5), 100*data["v_over_H0"])
    grid = np.linspace(s[0], CUT, 32769)
    eta = cumulative_simpson(np.exp(-loga(grid)), x=grid, initial=0.)
    return float(s[0]), (PchipInterpolator(eta, grid), loga, mass, eta[-1])


def source_covariance(q, a, mass, X, P, C):
    inertia = a**3*np.hypot(q/a, mass)
    number = (inertia*X + P/inertia - 1)/2
    corr = (inertia*X - P/inertia)/2 + 1j*C
    if min(number) < -1e-8 or np.min(number*(number+1) - abs(corr)**2) < -2e-8:
        raise ValueError("Native covariance violates the bosonic uncertainty bound")
    return number, corr


def direct_response(method, background, q, steps, order):
    # Exact-oscillator midpoint transfer and direct Gauss integration at ω>=q.
    # The per-q output is used without the uniform Simpson average in the
    # preceding helper: q here belongs to the native Gauss momentum rule.
    curves = method.base.modes(background, q, steps)
    output = {}
    labels = method.OUTPUTS
    for window, support in WINDOWS.items():
        integrals = method.integrated(curves, q, support, order)
        for band, result in integrals.items():
            per_q = result["per_q"]
            vacuum = per_q[:, labels.index("vacuum")]
            diagonal = per_q[:, labels.index("incoherent_excess")]
            real = per_q[:, labels.index("cross_real")]/np.sqrt(2)
            imag = per_q[:, labels.index("cross_imag")]/np.sqrt(2)
            output[window, band] = np.stack((vacuum, diagonal, real, imag), axis=1)
    return output


def fourier_response(start, background, q, samples, padding):
    # Independent proper-time propagation using the covariance's instantaneous
    # physical basis at the fork. Shift the Fourier grid to the exact threshold.
    clock, loga, mass, elapsed = background
    N = len(q)
    a0 = np.exp(loga(start))
    inertia = a0**3*np.hypot(q/a0, mass(start))
    chi0 = 1/np.sqrt(2*inertia)
    def rhs(s, state):
        a = np.exp(loga(s))
        freq2 = (q/a)**2+mass(s)**2
        return np.r_[state[N:]/a**3, -a**3*freq2*state[:N]]
    sol = solve_ivp(rhs, (start,CUT), np.r_[chi0,-1j*inertia*chi0], method="DOP853",
                    dense_output=True, rtol=2e-11, atol=1e-13, max_step=.005)
    if not sol.success:
        raise RuntimeError(sol.message)
    af = np.exp(loga(CUT))
    fine = np.linspace(0,elapsed,32769)
    inverse = PchipInterpolator(clock(fine), fine)
    output = {}
    nodes, gauss = np.polynomial.legendre.leggauss(96)
    for window, (lo,hi) in WINDOWS.items():
        left,right = map(float,inverse([lo,hi]))
        t = np.arange(samples)*((right-left)/samples)
        dt = (right-left)/samples
        proper = clock(left+t)
        wave = sol.sol(proper)[:N]*np.exp(loga(proper))*np.sin(np.pi*(proper-lo)/(hi-lo))**4
        nfft = samples*padding
        full, screened = [], []
        def products(pos, neg):
            return np.vstack((abs(neg)**2, abs(pos)**2+abs(neg)**2,
                              2*(pos*neg).real, 2*(pos*neg).imag))
        for k, g in zip(q,wave,strict=True):
            freq = k+np.arange(nfft//2)*2*np.pi/(nfft*dt)
            shift = np.exp(1j*k*t)
            plus = nfft*dt*np.fft.ifft(g*shift,n=nfft)[:nfft//2]
            minus = dt*np.fft.fft(g*shift.conj(),n=nfft)[:nfft//2]
            density = products(plus,minus)*((freq**2-k*k)/(8*np.pi**2))
            full.append(np.stack((simpson(density,x=freq), simpson(density*freq/af,x=freq))))
            w = k+(nodes+1)*(60-k)/2
            phase = np.exp(1j*w[:,None]*t)
            density = products(dt*(phase@g),dt*(phase.conj()@g))*(w*w-k*k)*gauss*(60-k)/(16*np.pi**2)
            screened.append(np.stack((np.sum(density,axis=1),np.sum(density*w/af,axis=1))))
        output[window,"full"] = np.stack(full,axis=-1)
        output[window,"screen"] = np.stack(screened,axis=-1)
    return output


def evaluate(response, number, corr, weight, af):
    result = {}
    for (window, band), matrix in response.items():
        vac, diag, real, imag = matrix.transpose(1,0,2)
        if np.min(vac) < 0 or np.min(diag-vac) < -1e-12:
            raise ValueError("Negative spectral branch")
        if np.max((real**2+imag**2)/4-vac*(diag-vac)) > 1e-8:
            raise ValueError("Emission cross term violates Cauchy–Schwarz")
        excess_modes = number*diag+corr.real*real-corr.imag*imag
        if np.min(vac+excess_modes) < -1e-10:
            raise ValueError("Negative pair yield")
        vacuum, extra = vac@weight, excess_modes@weight
        result.setdefault(window, {})[band] = {
            key: {"vacuum": float(vacuum[j]), "excess": float(extra[j]),
                  "total": float(vacuum[j]+extra[j]),
                  "total_density": float((vacuum[j]+extra[j])/af**3)}
            for j,key in enumerate(("number","energy"))
        }
    return result
