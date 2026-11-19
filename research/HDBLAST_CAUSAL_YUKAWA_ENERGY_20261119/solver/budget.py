"""Causal order-y² energy and work on the supplied expanding quiet wall."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicHermiteSpline


def read_native(native, grid, arm="quiet"):
    """Read one retained wall and the three alternative fork preparations."""
    if arm not in ("quiet", "plus", "minus", "far"):
        raise ValueError("Unknown native arm")
    native = Path(native)
    fork = json.loads((native / "FORK.json").read_bytes())[grid]
    with zipfile.ZipFile(native / "WALLS.zip") as archive:
        with archive.open(f"{grid}/retain/{arm}/wall.csv") as stream:
            wall = np.genfromtxt(stream, delimiter=",", names=True)
    s = wall["H0tau"]
    loga = CubicHermiteSpline(s, wall["A_b"] - fork["reference_A"],
                              wall["H_over_H0"], extrapolate=False)
    mass = CubicHermiteSpline(s, 100 * (wall["phi_b"] - .5),
                              100 * wall["v_over_H0"], extrapolate=False)
    start = float(s[0])
    q, weights = map(np.array, (fork["q"], fork["weight"]))
    nodes, dw = leggauss(len(q))
    np.testing.assert_allclose(q, 4 * (nodes + 1), rtol=1e-13)
    np.testing.assert_allclose(weights, 4 * dw * q**2 / (2 * np.pi**2), rtol=2e-12)
    np.testing.assert_allclose([start, np.exp(loga(start))],
                               [fork["H0tau"], fork["a"]], rtol=2e-13)
    states = {name: tuple(np.asarray(fork["states"][name][key]) for key in ("X", "C", "P"))
              for name in ("retain", "erase")}
    a = fork["a"]
    inertia = a**3 * np.hypot(q / a, mass(start))
    states["fork_vacuum"] = (1 / (2 * inertia), np.zeros_like(q), inertia / 2)
    return q, weights, states, (start, loga, mass)


# Small variables at each state and q. The static correction is evolved on its
# own and added to the bare solution, using linearity at this perturbative order.
SMALL = ("X", "C", "P", "X2", "C2", "P2", "Xct", "Cct", "Pct", "Dconf",
         "mass", "scale_parent", "switch", "scale_daughter", "Jparent", "Jdaughter",
         "mass_ct_response", "scale_ct_response", "switch_ct", "scale_ct", "Jctparent")
IDX = {name: i for i, name in enumerate(SMALL)}


def make_clock(start, window, end, steps):
    lo, hi = window
    middle = (lo + hi) / 2
    bounds = np.array([start, lo, middle, hi, end])
    if min(np.diff(bounds)) <= 0:
        raise ValueError("Require fork < switch on < switch off < cut")
    index = np.rint(steps * (bounds - start) / (end - start)).astype(int)
    if min(np.diff(index)) < 1:
        raise ValueError("Too few time steps for the support")
    grid = np.concatenate([np.linspace(left, right, n + 1)[:-1]
                           for left, right, n in zip(bounds[:-1], bounds[1:],
                                                     np.diff(index), strict=True)] + [bounds[-1:]])
    return grid, index


def run_budget(q, weights, initial, background, cutoff, order=64, steps=1024,
               window=(1.35, 2.05), end=2.10, history=False):
    """Integrate bare and statically subtracted budgets; weights may select one q."""
    start, loga, mass = background
    q, weights = np.asarray(q), np.asarray(weights)
    if not (0 < min(q) <= max(q) < cutoff) or q.shape != weights.shape:
        raise ValueError("Require positive momenta below the comoving bath cutoff")
    at, sample = make_clock(start, window, end, steps)
    time = np.empty(2 * steps + 1)
    time[::2] = at
    time[1::2] = (at[:-1] + at[1:]) / 2
    a = np.exp(loga(time))
    h = loga(time, 1)
    m = mass(time)
    dm = mass(time, 1)
    if not np.all(np.isfinite(a + h + m + dm)):
        raise ValueError("Native history does not cover the requested evolution")
    lo, hi = window
    angle = np.pi * (time - lo) / (hi - lo)
    active = (time > lo) & (time < hi)
    f = np.where(active, np.sin(angle)**4, 0.)
    df = np.where(active, 4 * np.pi / (hi - lo) * np.sin(angle)**3 * np.cos(angle), 0.)
    inv_inertia = a**-3
    stiffness = a[:, None] * q**2 + (a**3 * m**2)[:, None]
    scale_coeff = .5 * h[:, None] * (a[:, None] * q**2 + (3 * a**3 * m**2)[:, None])
    mass_coeff = a**3 * m * dm

    node, quadrature = leggauss(order)
    omega = q[:, None] + (cutoff - q[:, None]) * (node + 1) / 2
    domega = (cutoff - q[:, None]) * quadrature / 2
    omega2 = omega**2
    g2 = omega * (omega2 - q[:, None]**2) * domega / (4 * np.pi**2)
    ct = np.sum(g2 / omega2, axis=-1)
    labels = tuple(initial)
    small = np.zeros((len(SMALL), len(labels), len(q)))
    small[:3] = np.asarray([initial[k] for k in labels]).transpose(1, 0, 2)
    # Each mixed covariance is divided by y*g_j. Their unperturbed values vanish.
    mix = np.zeros((4, len(labels), len(q), order))

    def weighted_sum(v):
        return np.einsum("sqj,qj->sq", v, g2, optimize=False)

    def rhs(i, z, t):
        X, C, P, X2, C2, P2, Xc, Cc, Pc, D = z[:10]
        Qx, Qp, Px, Pp = t
        A, B, F, ai = inv_inertia[i], stiffness[i], f[i], a[i]
        V = ai * F**2 * ct
        mx, mp, px = map(weighted_sum, (Qx, Qp, Px))
        dz = np.empty_like(z)
        dz[:3] = (2 * A * C, A * P - B * X, -2 * B * C)
        dz[3:6] = (2 * A * C2, A * P2 - B * X2 - F * mx, -2 * (B * C2 + F * px))
        dz[6:9] = (2 * A * Cc, A * Pc - B * Xc - V * X, -2 * (B * Cc + V * C))
        dz[9] = -F * mp
        # Explicit parameter derivatives of the Hamiltonian at fixed canonical
        # variables. No work is inferred from an endpoint difference.
        dz[10] = mass_coeff[i] * X2
        dz[11] = scale_coeff[i] * X2 - 1.5 * h[i] * A * P2
        dz[12] = df[i] * mx
        dz[13] = -h[i] * D / ai
        dz[14] = -F * A * px
        dz[15] = -F * mp / ai
        dz[16] = mass_coeff[i] * Xc
        dz[17] = scale_coeff[i] * Xc - 1.5 * h[i] * A * Pc
        dz[18] = ai * F * df[i] * ct * X
        dz[19] = .5 * h[i] * V * X
        dz[20] = -V * A * C
        dt = np.empty_like(t)
        dt[0] = A * Px + Qp / ai
        dt[1] = A * Pp - omega2 * Qx / ai - F * X[..., None]
        dt[2] = -B[:, None] * Qx + Pp / ai - F / (2 * omega)
        dt[3] = -B[:, None] * Qp - omega2 * Px / ai - F * C[..., None]
        return dz, dt

    def ledger(i, z, t):
        A, B, ai, F = inv_inertia[i], stiffness[i], a[i], f[i]
        v = {k: z[j] for j, k in enumerate(SMALL)}
        e_d = v["Dconf"] / ai
        e_int = F * weighted_sum(t[0])
        rows = {}
        for scheme, use in (("bare", 0), ("static", 1)):
            X2 = v["X2"] + use * v["Xct"]
            P2 = v["P2"] + use * v["Pct"]
            e_p = .5 * (A * P2 + B * X2)
            e_ct = .5 * use * ai * F**2 * ct * v["X"]
            energy = dict(parent=e_p, daughter=e_d, interaction=e_int, counterterm=e_ct)
            work = dict(mass=v["mass"] + use * v["mass_ct_response"],
                        scale_parent=v["scale_parent"] + use * v["scale_ct_response"],
                        scale_daughter=v["scale_daughter"], scale_counterterm=use * v["scale_ct"],
                        switch=v["switch"], switch_counterterm=use * v["switch_ct"])
            transfer = dict(int_to_parent=v["Jparent"], int_to_daughter=v["Jdaughter"],
                            ct_to_parent=use * v["Jctparent"])
            residual = dict(total=sum(energy.values()) - sum(work.values()),
                            parent=e_p-work["mass"]-work["scale_parent"]
                                   -transfer["int_to_parent"]-transfer["ct_to_parent"],
                            daughter=e_d-work["scale_daughter"]-transfer["int_to_daughter"],
                            interaction=e_int-work["switch"]+transfer["int_to_parent"]
                                        +transfer["int_to_daughter"],
                            counterterm=e_ct-work["switch_counterterm"]-work["scale_counterterm"]
                                        +transfer["ct_to_parent"])
            each = {key: {field: val @ weights for field, val in section.items()}
                    for key, section in (("energy", energy), ("work", work),
                                         ("transfer", transfer), ("residual", residual))}
            rows[scheme] = {label: {key: {field: float(val[j]) for field, val in section.items()}
                                    for key, section in each.items()} for j, label in enumerate(labels)}
        return rows

    names = ("fork", "switch_on", "peak", "switch_off", "cut")
    saves = dict(zip(sample, names, strict=True))
    history_at = set(np.rint(np.linspace(0, steps, 65)).astype(int)) | set(saves)
    result = {k: {scheme: {"checkpoints": {}, "max_abs_residual": {}}
                  for scheme in ("bare", "static")} for k in labels}

    def save(step):
        rows = ledger(step * 2, small, mix)
        for scheme, states in rows.items():
            for label, row in states.items():
                rec = result[label][scheme]
                for k, x in row["residual"].items():
                    rec["max_abs_residual"][k] = max(rec["max_abs_residual"].get(k, 0.), abs(x))
                if step in saves:
                    rec["checkpoints"][saves[step]] = dict(s=float(at[step]), a=float(a[step*2]), **row)
                if history:
                    rec.setdefault("trace", []).append(dict(s=float(at[step]), **row))

    save(0)
    for k, delta in enumerate(np.diff(at)):
        i = k * 2
        s1, t1 = rhs(i, small, mix)
        s2, t2 = rhs(i+1, small + delta/2*s1, mix + delta/2*t1)
        s3, t3 = rhs(i+1, small + delta/2*s2, mix + delta/2*t2)
        s4, t4 = rhs(i+2, small + delta*s3, mix + delta*t3)
        small += delta/6 * (s1 + 2*s2 + 2*s3 + s4)
        mix += delta/6 * (t1 + 2*t2 + 2*t3 + t4)
        if k+1 in history_at:
            save(k+1)
    return result


def all_runs(native):
    native = Path(native)
    result = dict(schema=1, cutoff_frame="comoving", units="y^2 H0^4 per comoving volume",
                  inputs={p: hashlib.sha256((native/p).read_bytes()).hexdigest()
                          for p in ("FORK.json", "WALLS.zip")}, runs={})
    for grid in ("fine", "coarse"):
        inputs = read_native(native, grid)
        for window, bounds in (("wide", (1.35, 2.05)), ("narrow", (1.45, 1.95))):
            cutoffs = (40, 60, 80, 120) if window == "wide" and grid == "fine" else (60, 120)
            for cutoff in cutoffs:
                key = f"{grid}_{window}_{cutoff}"
                result["runs"][key] = run_budget(*inputs, cutoff, window=bounds,
                                                  history=(key=="fine_wide_60"))
    # Separate Δs and frequency-order refinements at the most oscillatory cutoff.
    inputs = read_native(native, "fine")
    for order, steps in ((64, 2048), (96, 1024), (96, 2048)):
        key = f"fine_wide_120_N{order}_T{steps}"
        result["runs"][key] = run_budget(*inputs, 120, order, steps)
    q, weights, initial, bg = inputs
    retain = {"retain": initial["retain"]}
    for cutoff, order, steps in ((240, 128, 2048), (240, 160, 4096),
                                 (480, 192, 4096), (480, 240, 8192)):
        key = f"fine_wide_{cutoff}_N{order}_T{steps}"
        result["runs"][key] = run_budget(q, weights, retain, bg, cutoff, order, steps)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native", type=Path, help="Prior capsule's native directory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = all_runs(args.native)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
