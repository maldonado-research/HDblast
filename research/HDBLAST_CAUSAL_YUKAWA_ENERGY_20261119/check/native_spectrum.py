"""Independent proper-time outgoing daughter-energy calculation.

No preceding response code is imported. Read the public native wall and fork
archives directly, solve the parent mode and conformal clock with DOP853, then
integrate the frequency amplitudes by Simpson's rule on a proper-time grid.
"""
import argparse
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline


def read_native(native, grid="fine", geometry="retain", arm="quiet"):
    native = Path(native)
    fork = json.loads((native / "FORK.json").read_text())[grid]
    with zipfile.ZipFile(native / "WALLS.zip") as archive:
        raw = archive.read(f"{grid}/{geometry}/{arm}/wall.csv")
    wall = np.genfromtxt(io.BytesIO(raw), names=True, delimiter=",")
    time = wall["H0tau"]
    loga = CubicHermiteSpline(time, wall["A_b"] - fork["reference_A"], wall["H_over_H0"])
    mass = CubicHermiteSpline(time, 100*(wall["phi_b"]-.5), 100*wall["v_over_H0"])
    q = np.asarray(fork["q"])
    weight = np.asarray(fork["weight"])
    initial = {name: tuple(np.asarray(row[key]) for key in ("X", "C", "P"))
               for name, row in fork["states"].items() if name in ("retain", "erase")}
    inertia = fork["a"]**3*np.hypot(q/fork["a"], 100*(fork["phi_b"]-.5))
    initial["fork_vacuum"] = (1/(2*inertia), np.zeros_like(q), inertia/2)
    if not np.isclose(np.exp(loga(time[0])), fork["a"], rtol=2e-10):
        raise ValueError("The wall and initial covariance use different scale-factor conventions")
    if time[-1] < 2.10:
        raise ValueError("Wall does not reach readout")
    return q, weight, initial, (float(time[0]), loga, mass)


def spectrum(q, weights, initial, background, cutoff=60., order=96,
             samples=4096, window=(1.35, 2.05), end=2.10):
    start, loga, mass = background
    q, weights = np.asarray(q), np.asarray(weights)
    if not start < window[0] < window[1] < end:
        raise ValueError("Switch must vanish before preparation and readout")
    if cutoff <= max(q):
        raise ValueError("Cutoff below a parent momentum")
    a0 = np.exp(loga(start))
    inertia = a0**3*np.hypot(q/a0, mass(start))
    u0 = 1/np.sqrt(2*inertia)
    nq = len(q)
    def rhs(t, state):
        a = np.exp(loga(t))
        return np.r_[state[nq:2*nq]/a**3,
                     -a**3*((q/a)**2+mass(t)**2)*state[:nq], 1/a]
    solution = solve_ivp(rhs, (start, end), np.r_[u0, -1j*inertia*u0, 0j],
                         method="DOP853", rtol=2e-12, atol=2e-14, dense_output=True,
                         max_step=.002)
    if not solution.success:
        raise RuntimeError(solution.message)
    t = np.linspace(*window, samples+1)
    curves = solution.sol(t)
    switch = np.sin(np.pi*(t-window[0])/(window[1]-window[0]))**4
    nodes, quadrature = np.polynomial.legendre.leggauss(order)
    af = float(np.exp(loga(end)))
    per_q = {name: np.zeros((2, nq)) for name in initial}
    decomposition = {name: np.zeros((3, 2, nq)) for name in initial}
    for j, momentum in enumerate(q):
        omega = momentum+(nodes+1)*(cutoff-momentum)/2
        dw = quadrature*(cutoff-momentum)/2
        phase = np.exp(1j*omega[:, None]*curves[-1].real)
        wave = switch*curves[j]
        positive = simpson(phase*wave, x=t, axis=1)
        negative = simpson(phase.conj()*wave, x=t, axis=1)
        density = (omega**2-momentum**2)/(8*np.pi**2)*dw
        kernels = np.stack([np.ones_like(omega), omega/af]) * density
        for name, (X, C, P) in initial.items():
            number = (inertia[j]*X[j]+P[j]/inertia[j]-1)/2
            corr = (inertia[j]*X[j]-P[j]/inertia[j])/2+1j*C[j]
            pieces = np.array([abs(negative)**2,
                               number*(abs(positive)**2+abs(negative)**2),
                               2*np.real(corr*positive*negative)])
            decomposition[name][:, :, j] = pieces@kernels.T
            per_q[name][:, j] = pieces.sum(axis=0)@kernels.T
    output = {}
    for name, values in per_q.items():
        total = values@weights
        parts = decomposition[name]@weights
        output[name] = dict(number=float(total[0]), energy=float(total[1]),
                            energy_density=float(total[1]/af**3),
                            vacuum_energy=float(parts[0, 1]), diagonal_energy=float(parts[1, 1]),
                            coherent_energy=float(parts[2, 1]),
                            per_q_energy=values[1].tolist(), per_q_number=values[0].tolist())
    return dict(cutoff=float(cutoff), frequency_order=order, time_intervals=samples,
                proper_window=list(window), proper_readout=end, a_final=af, states=output)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--cutoffs", nargs="+", type=float, default=[40, 60, 80, 120])
    ap.add_argument("--grids", nargs="+", default=["coarse", "fine"])
    ap.add_argument("--order", type=int, default=96)
    ap.add_argument("--samples", type=int, default=4096)
    args = ap.parse_args()
    cutoff_key = lambda cut: str(int(cut)) if float(cut).is_integer() else str(cut)
    result = dict(method="independent DOP853 parent and proper-time Simpson amplitudes", grids={}, narrow_grids={})
    for grid in args.grids:
        inputs = read_native(args.native, grid)
        result["grids"][grid] = {cutoff_key(cut): spectrum(*inputs, cutoff=cut, order=args.order,
                                                        samples=args.samples) for cut in args.cutoffs}
        result["narrow_grids"][grid] = {cutoff_key(cut): spectrum(*inputs, cutoff=cut, order=args.order,
                                                              samples=args.samples, window=(1.45, 1.95))
                                         for cut in args.cutoffs if cut in (60., 120.)}
    args.output.parent.mkdir(exist_ok=True, parents=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
