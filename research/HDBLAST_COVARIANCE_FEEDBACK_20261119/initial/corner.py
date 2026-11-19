"""Check the time derivative of the shell junctions in a native covariance fork."""
import argparse
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline


HERE = Path(__file__).resolve().parent
WALLS = Path("research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip")
WALL_SHA = "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"
G = 100.
RB = 7.838285538073243
B = 6.49e-6


def checked(path, digest):
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != digest:
        raise ValueError(f"Unexpected SHA-256 for {path}")
    return payload


def quadrature(count):
    x, w = np.polynomial.legendre.leggauss(count)
    q = 4*(1+x)
    return q, 4*w*q*q/(2*np.pi**2)


def background(archive, grid, ref):
    def records(arm):
        return np.genfromtxt(archive.open(f"{grid}/y0_b1/{arm}.csv"), delimiter=",", names=True)
    prefix, quiet = records("pre"), records("quiet")
    fields = np.concatenate((prefix, quiet))
    if np.any(np.diff(fields["H0tau"]) <= 0):
        raise ValueError("Original wall clock does not increase")
    loga = CubicHermiteSpline(fields["H0tau"], fields["A_b"] - ref, fields["H_over_H0"])
    phi = CubicHermiteSpline(fields["H0tau"], fields["phi_b"], fields["v_over_H0"])
    return loga, phi, prefix[0]["H0tau"], quiet[0]["H0tau"]


def covariance_rhs(clock, y, loga, phi, q):
    x, p, c = y.reshape(3, -1)
    a = np.exp(loga(clock))
    e2 = (q/a)**2 + (G*(phi(clock)-.5))**2
    return np.concatenate((2*c/a**3, -2*a**3*e2*c, p/a**3-a**3*e2*x))


def reconstruct(funs, source, cut, q):
    loga, phi = funs
    a = np.exp(loga(source))
    e = np.hypot(q/a, G*(phi(source)-.5))
    start = np.concatenate((.5/(a**3*e), .5*a**3*e, np.zeros_like(q)))
    sol = solve_ivp(covariance_rhs, (source, cut), start, args=(loga,phi,q),
                    rtol=2e-11, atol=2e-13, method="DOP853", max_step=.0008)
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol.y[:,-1].reshape(3, -1)


def smooth_cut(t, w, c):
    basis = np.stack((np.ones_like(t), t), axis=1)
    moments = w*c
    b0, b1 = np.linalg.solve(basis.T @ (moments[:, None]*basis), basis.T @ (moments*t*t))
    test = np.array([0., 1., np.clip(b1/2., 0., 1.)])
    norm = float(max(abs(test**2-b1*test-b0)))
    return .9*c*(t*t-b1*t-b0)/norm, [float(b0), float(b1), norm]


def cuts(funs, fork_time, q, w, a, phi, delta_c):
    # Propagate the difference of covariance matrices. The two preparations
    # share X and P at the fork, so their initial differences vanish modewise.
    loga_s, phi_s = funs
    times = np.array([.00075,.0015,.003,.006,.012,.024,.05,.1])
    elapsed = times+fork_time
    initial = np.concatenate((np.zeros_like(q), np.zeros_like(q), delta_c))
    sol = solve_ivp(covariance_rhs, (fork_time,elapsed[-1]), initial, args=(loga_s,phi_s,q),
                    t_eval=elapsed, rtol=2e-12, atol=2e-15, method="DOP853", max_step=.00025)
    if not sol.success:
        raise ArithmeticError(sol.message)
    m2 = w @ (q**4*delta_c)
    p3, j5 = 16*m2/(3*a**7), 32*G**2*(phi-.5)*m2/a**7
    curve = []
    for dt, s, y in zip(times, elapsed, sol.y.T):
        a_s, delta = np.exp(loga_s(s)), phi_s(s)-.5
        x, p, _ = y.reshape(3,-1)
        pressure = .5*w @ (p/a_s**6-((G*delta)**2+(q/a_s)**2/3)*x)
        force = G**2*delta*(w @ x)
        curve.append(dict(elapsed=float(dt),pressure=float(pressure),force=float(force),
                          pressure_over_cubic=float(pressure/(p3*dt**3/6)),
                          force_over_quintic=float(force/(j5*dt**5/120))))
    return dict(pressure_third=float(p3),force_fifth=float(j5),curve=curve)


def run(row, geom):
    loga, phi_s, source, cut = geom
    x, p, c = np.array(row["covariance"])
    q, w = quadrature(len(x))
    a, phi = float(np.exp(loga(cut))), float(phi_s(cut))
    if abs(cut-row["s"])>1e-12 or abs(a-row["a"])>1e-12 or abs(phi-row["phi"])>1e-12:
        raise ValueError("Wall and fork initial values disagree")
    e2 = (q/a)**2+(G*(phi-.5))**2
    cap = np.sqrt(x*p-.25)
    if min(c)<=0 or min(x*p-c*c)<.25-2e-10:
        raise ValueError("Expected a positive, initially pure covariance orbit")
    # With exact equality of ∫w C, downward and upward changes have equal
    # integrated magnitude. Positive room above the native C controls both.
    slack = np.maximum(0.,cap-c)
    residual_bound = float((w@slack)/(w@c))
    read = {}
    for name, diff in (("retain",np.zeros_like(c)),("reflect",-2*c),("dephase",-c)):
        m0, m1 = w@diff, w@(q*q*diff)
        p_prime = -2/a**3*((G*(phi-.5))**2*m0+2*m1/(3*a*a))
        j_prime = 2*G**2*(phi-.5)*m0/a**3
        factor = -.5*RB*row["lapse"]**2*B/RB
        read[name]=dict(pressure_rate_change=float(p_prime),force_rate_change=float(j_prime),
                        B_TZ_change=float(factor*p_prime),phi_TZ_change=float(factor*j_prime))
    smooth, coef = smooth_cut((q/8)**2,w,c)
    new_delta = 2*smooth
    m0, m1 = w@new_delta, w@(q*q*new_delta)
    if max(abs(m0),abs(m1))>1e-12:
        raise ValueError("Moment matching failed")
    # The same wall supplies a quadrature-independent free-mode replay,
    # at the native nodes as well as at 257 uniformly spaced momenta.
    at_native = reconstruct((loga,phi_s),source,cut,q)
    probes = np.linspace(0.,8.,257)
    at_probes = reconstruct((loga,phi_s),source,cut,probes)[2]
    return dict(C_min=float(c.min()), C_max=float(c.max()),
                max_pure_determinant_error=float(max(abs(x*p-c*c-.25))),
                max_relative_same_force_derivative_perturbation=residual_bound,
                native_integral_C=float(w@c), native_integral_q2_C=float(w@(q*q*c)),
                independent_max_C_error_native_nodes=float(max(abs(c-at_native[2]))),
                independent_257_min_C=float(min(at_probes)), independent_minimizing_q=float(probes[at_probes.argmin()]),
                junction_derivatives=read,
                smooth_control=dict(polynomial_coefficients=coef, delta_C_moment0=float(m0),
                                    delta_C_moment1=float(m1),delta_C_moment2=float(w@(q**4*new_delta)),
                                    min_determinant=float(min(x*p-smooth**2)),
                                    **cuts((loga,phi_s),cut,q,w,a,phi,new_delta)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo",type=Path)
    opts = parser.parse_args()
    forks = json.loads((HERE/"FORK.json").read_bytes())
    with ZipFile(io.BytesIO(checked(opts.repo/WALLS,WALL_SHA))) as archive:
        out = {grid:run(row,background(archive,grid,row["log_a_reference"])) for grid,row in forks.items()}
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
