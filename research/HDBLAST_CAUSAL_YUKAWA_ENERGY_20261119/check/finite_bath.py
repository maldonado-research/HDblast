"""Independent symplectic finite-bath energy and weak-y extrapolation check.

The parent and a quadrature-discretized oscillator bath form one exact quadratic
Hamiltonian. A midpoint kick-drift-kick map preserves the canonical symplectic
form at every step. Propagate its full transfer matrix at finite y and propagate
the first three Taylor coefficients of the same map independently. No solver
hierarchy is imported. The map error is checked by halving its proper-time step.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

from native_spectrum import read_native, spectrum


def switch(t, window):
    lo, hi = window
    if not lo < t < hi:
        return 0., 0.
    phase = np.pi*(t-lo)/(hi-lo)
    s, c = np.sin(phase), np.cos(phase)
    return s**4, 4*np.pi/(hi-lo)*s**3*c


def symplectic_run(q, initial, background, cutoff=60., order=24, steps=8192,
                   window=(1.35, 2.05), end=2.10, static=False,
                   couplings=(0., .05, .025)):
    start, loga, mass = background
    z, quadrature = np.polynomial.legendre.leggauss(order)
    omega = q+(z+1)*(cutoff-q)/2
    dw = quadrature*(cutoff-q)/2
    g = np.sqrt(2*omega*(omega**2-q*q)/(8*np.pi**2)*dw)
    counterterm = float(np.sum(g*g/(omega*omega))) if static else 0.
    n = order+1
    finite = np.repeat(np.eye(2*n)[None], len(couplings), axis=0)
    series = np.zeros((3, 2*n, 2*n))
    series[0] = np.eye(2*n)
    X, C, P = map(float, initial)
    covariance_diagonal = np.r_[X, 1/(2*omega), P, omega/2]
    def cov(L, R):
        # L C_initial R^T, including the parent's X-P correlation.
        return (L*covariance_diagonal)@R.T+C*(np.outer(L[:, 0], R[:, n])+
                                                      np.outer(L[:, n], R[:, 0]))
    def coefficient_covariance():
        zero = cov(series[0], series[0])
        one = cov(series[0], series[1])+cov(series[1], series[0])
        two = cov(series[0], series[2])+cov(series[2], series[0])+cov(series[1], series[1])
        return zero, one, two
    def diagonal_cov(L, R):
        return np.sum(L*R*covariance_diagonal, axis=-1)+C*(L[:, 0]*R[:, n]+L[:, n]*R[:, 0])
    def series_observables(t):
        a = np.exp(loga(t)); H = loga(t, 1)
        m = mass(t); mdot = mass(t, 1)
        f, fdot = switch(t, window)
        xx0 = diagonal_cov(series[0, :n], series[0, :n])
        xx2 = (2*diagonal_cov(series[0, :n], series[2, :n])+
               diagonal_cov(series[1, :n], series[1, :n]))
        pp2 = (2*diagonal_cov(series[0, n:], series[2, n:])+
               diagonal_cov(series[1, n:], series[1, n:]))
        cross1 = cov(series[0, :1], series[1, 1:n])+cov(series[1, :1], series[0, 1:n])
        G = float(cross1[0]@g)
        parent = .5*(pp2[0]/a**3+a**3*((q/a)**2+m*m)*xx2[0])
        bath = .5*np.sum(pp2[1:]+omega*omega*xx2[1:])/a
        interaction = f*G
        ct = .5*counterterm*a*f*f*xx0[0]
        work_parent = .5*(-3*H*pp2[0]/a**3+
                          a**3*(H*((q/a)**2+3*m*m)+2*m*mdot)*xx2[0])
        work_bath = -H*bath
        work_switch = fdot*G
        work_ct = .5*counterterm*a*(H*f*f+2*f*fdot)*xx0[0]
        return np.array([parent, bath, interaction, ct,
                         work_parent, work_bath, work_switch, work_ct])
    dt = (end-start)/steps
    times = np.linspace(start, end, steps+1)
    trace = np.empty((steps+1, 8))
    trace[0] = series_observables(start)
    couplings = np.asarray(couplings)[:, None, None]
    for j in range(steps):
        t = (times[j]+times[j+1])/2
        a = np.exp(loga(t)); m = mass(t); f = switch(t, window)[0]
        kinetic = np.r_[a**-3, np.full(order, 1/a)]
        potential = np.r_[a**3*((q/a)**2+m*m), omega*omega/a]
        def cross(coordinates):
            output = np.empty_like(coordinates)
            output[..., 0, :] = f*np.einsum("i,...ij->...j", g, coordinates[..., 1:, :])
            output[..., 1:, :] = f*g[:, None]*coordinates[..., :1, :]
            return output
        def kick(duration):
            x = finite[:, :n]
            finite[:, n:] -= duration*(potential[None, :, None]*x+couplings*cross(x))
            finite[:, n] -= duration*couplings[:, 0, 0, None]**2*counterterm*a*f*f*x[:, 0]
            # Polynomial multiplication includes A_1 M_1 in the y^2 term.
            sx = series[:, :n]
            force = potential[None, :, None]*sx
            force[1:] += cross(sx[:-1])
            force[2, 0] += counterterm*a*f*f*sx[0, 0]
            series[:, n:] -= duration*force
        kick(dt/2)
        finite[:, :n] += dt*kinetic[None, :, None]*finite[:, n:]
        series[:, :n] += dt*kinetic[None, :, None]*series[:, n:]
        kick(dt/2)
        trace[j+1] = series_observables(times[j+1])
    exact = []
    a = np.exp(loga(end)); m = mass(end)
    for coupling, transfer in zip(couplings[:, 0, 0], finite, strict=True):
        xx = diagonal_cov(transfer[:n], transfer[:n])
        pp = diagonal_cov(transfer[n:], transfer[n:])
        exact.append(dict(y=float(coupling),
                          parent=float(.5*(pp[0]/a**3+a**3*((q/a)**2+m*m)*xx[0])),
                          bath=float(.5*np.sum(pp[1:]+omega*omega*xx[1:])/a)))
    exact_coefficients = []
    for row in exact[1:]:
        exact_coefficients.append({key: (row[key]-exact[0][key])/row["y"]**2 for key in ("parent", "bath")})
    ratio = (exact[1]["y"]/exact[2]["y"])**2
    extrapolated = {key: (ratio*exact_coefficients[1][key]-exact_coefficients[0][key])/(ratio-1)
                    for key in ("parent", "bath")}
    works = simpson(trace[:, 4:], x=times, axis=0)
    initial_energy = float(trace[0, :4].sum())
    final_energy = float(trace[-1, :4].sum())
    endpoint = dict(zip(("parent", "bath", "interaction", "counterterm"), map(float, trace[-1, :4])))
    zero, one, two = coefficient_covariance()
    canonical = np.block([[np.zeros((n,n)), np.eye(n)], [-np.eye(n), np.zeros((n,n))]])
    symplectic = max(float(np.max(abs(T@canonical@T.T-canonical))) for T in finite)
    commutator2 = series[0]@canonical@series[2].T+series[2]@canonical@series[0].T+series[1]@canonical@series[1].T
    return dict(parent_momentum=float(q), cutoff=cutoff, frequency_order=order,
                proper_time_steps=steps, static_counterterm=static,
                finite_coupling_endpoints=exact, finite_y_coefficients=exact_coefficients,
                weak_y_extrapolated=extrapolated, order_y2_endpoint=endpoint,
                extrapolation_minus_series={key: extrapolated[key]-endpoint[key] for key in extrapolated},
                order_y2_work=dict(zip(("parent_background", "bath_redshift", "switch", "counterterm"), map(float, works))),
                energy_balance_residual=final_energy-initial_energy-float(works.sum()),
                max_symplectic_residual=symplectic,
                order_y2_commutator_residual=float(np.max(abs(commutator2))),
                initial_parent_uncertainty_determinant=X*P-C*C,
                min_baseline_covariance_eigenvalue_at_cut=float(min(np.linalg.eigvalsh(zero))))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--steps", nargs="+", type=int, default=[4096, 8192])
    ap.add_argument("--cutoff", type=float, default=60.)
    ap.add_argument("--order", type=int, default=24)
    ap.add_argument("--index", type=int, default=11)
    args = ap.parse_args()
    q, weights, initial, background = read_native(args.native)
    j = args.index
    inputs = {name: tuple(np.array([part[j]]) for part in triplet) for name, triplet in initial.items()}
    spectral = spectrum(q[j:j+1], np.ones(1), inputs, background,
                        cutoff=args.cutoff, order=args.order, samples=8192)
    runs = {}
    for name, triplet in initial.items():
        runs[name] = {}
        for static in (False, True):
            label = "static" if static else "bare"
            runs[name][label] = []
            for steps in args.steps:
                result = symplectic_run(q[j], tuple(x[j] for x in triplet), background,
                                        cutoff=args.cutoff, order=args.order, steps=steps, static=static)
                reference = spectral["states"][name]["energy"]
                result["spectrum_bath_energy"] = reference
                result["series_minus_spectrum"] = result["order_y2_endpoint"]["bath"]-reference
                runs[name][label].append(result)
    output = dict(method="exact finite-y symplectic transfer and independent Taylor expansion", grid="fine",
                  native_momentum_index=j, runs=runs)
    args.output.parent.mkdir(exist_ok=True, parents=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
