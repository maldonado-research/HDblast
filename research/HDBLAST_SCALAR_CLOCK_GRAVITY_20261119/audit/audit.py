"""Independent integrated scalar-flux and brane-expansion tests on bulk replay samples."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
CUT = 2.10


def integral(clock, values, stop, stride=1):
    x, y = clock[::stride], np.asarray(values)[::stride]
    if not x[0] < stop <= x[-1]:
        raise ValueError("The quadrature grid does not cover the readout")
    spline = CubicSpline(x, y)
    return float(spline.integrate(x[0], stop))


def at(clock, values, stop):
    if not clock[0] <= stop <= clock[-1]:
        raise ValueError("The trajectory does not cover the readout")
    return float(CubicSpline(clock, values)(stop))


def interior_gradient(z, fields, degree):
    """Extrapolate a derivative from strictly interior, possibly nonuniform, samples."""
    z = z[:degree + 1]
    if np.any(z >= 0) or np.any(np.diff(z) >= 0):
        raise ValueError("Interior nodes must be negative and ordered toward the bulk")
    length = abs(z[0])
    powers = np.vander(z / length, N=degree + 1, increasing=True).T
    target = np.eye(degree + 1)[1] / length
    weights = np.linalg.solve(powers, target)
    # Subtract the nearest node to avoid unnecessarily amplifying the offset.
    return (fields[:, :degree + 1] - fields[:, :1]) @ weights


def measure(d):
    s, t = d["s"], d["T"]
    if not (np.all(np.diff(s) > 0) and np.all(np.diff(t) > 0)):
        raise ValueError("Nonmonotone clock")
    a = np.exp(d["A"] - float(d["Aref"]))
    b, rb, rb2 = float(d["b"]), float(d["rb"]), float(d["rb"])**2
    mult = b / rb
    rho, p, j = (d[k] * mult for k in ("rho", "p", "j"))
    sig, sig1 = d["sigma"], d["sigma1"]
    phi, h, v = d["phi"], d["h"], d["v"]
    n = np.exp(-d["B"])
    vol = a**3
    tc = at(s, t, CUT)

    # Tetrad contraction of the scalar stress from both sides, with proper shell
    # volume. The orientation is independently fixed by the inward coordinate.
    tension = vol * sig1 * d["phit"]
    target = vol * j * d["phit"]
    individual = {"matter": integral(t, target, tc),
                  "work_rk": mult * (at(s, d["work"], CUT) - d["work"][0]),
                  "tension": integral(t, tension, tc)}
    for degree in (2, 4):
        grad = interior_gradient(d["z_inside"], d["inside_phi"], degree)
        rate = -2 * vol * n * d["phit"] * grad
        individual["bulk_" + str(degree)] = integral(t, rate, tc)
        for stride in (1, 2, 4):
            label = f"interior_{degree}" + (f"_stride_{stride}" if stride > 1 else "")
            individual[label] = integral(t, rate - tension, tc, stride)
    individual["matter_energy_change"] = mult * (at(s, vol*d["rho"], CUT) - vol[0]*d["rho"][0])
    individual["pressure_work"] = mult * integral(s, 3*vol*h*d["p"], CUT)
    individual["matter_balance_error"] = individual["matter_energy_change"] + individual["pressure_work"] - individual["matter"]
    individual["interior_error"] = individual["interior_4"] - individual["matter"]
    individual["quadratic_error"] = individual["interior_2"] - individual["matter"]

    # Integrate the Gauss/acceleration equation without constructing the Weyl
    # term from the Hubble rate. U is evaluated in superpotential form.
    W = 1. - phi + phi**3/3
    U = .5*(phi**2-1)**2 - (2/3)*W**2
    terms = {
        "vacuum": rb2*(sig**2/18 - sig1**2/24 + U/3),
        "scalar_kinetic": -v**2/6,
        "linear_stress_trace": rb2*sig*(rho - 3*p)/36,
        "quadratic_stress": -rb2*rho*(rho + 3*p)/36,
        "scalar_current": -rb2*j*(2*sig1 + j)/24,
    }
    a_final = at(s, a, CUT)
    weighted = lambda values: integral(s, a**2*values, CUT) / a_final**2
    trace = {key: weighted(term) for key, term in terms.items()}
    trace["initial"] = float(h[0]*(a[0]/a_final)**2)
    trace["predicted"] = sum(trace.values())
    # At fixed remaining matter the pure radiation (x,x/3) cancels from
    # sigma*(rho-3*p). Keep its exact part of the quadratic term separately.
    radiation = mult * d["rad"]
    rest_rho, rest_pressure = rho - radiation, p - radiation/3
    direct_rad = -rb2*((rest_rho + rest_pressure)*radiation/12 + radiation**2/18)
    trace["direct_radiation_part_of_quadratic_stress"] = weighted(direct_rad)
    trace["measured"] = at(s, h, CUT)
    trace["error"] = trace["measured"] - trace["predicted"]
    ga, gb, gf = (n[:, None]*d["actual_grad"]).T
    ja, jb, jf = (n[:, None]*d["junction"]).T
    grad_defect = (ga - ja)*(ga + ja + gb) + ja*(gb-jb) - (gf-jf)*(gf+jf)/6
    trace["stencil_correction"] = weighted(rb2*grad_defect)
    trace["hamiltonian_correction"] = weighted(-rb2*n**2*d["Hraw"]/6)
    trace["remainder"] = trace["error"] - trace["stencil_correction"] - trace["hamiltonian_correction"]
    phys_junctions = np.stack(((sig + rho)/6, (sig - 2*rho - 3*p)/6, -(sig1+j)/2), axis=1)
    maximum_junction_error = np.max(abs(n[:, None]*d["junction"]-phys_junctions))
    if maximum_junction_error > 2e-13:
        raise ValueError("Unexpected physical signs or matter-unit normalization at junctions")
    if abs(at(t, s, tc)-CUT) > 1e-11:
        raise ValueError("Proper/conformal clock inversion")
    if abs(trace["remainder"]) > 1e-4 or abs(individual["matter"]-individual["work_rk"]) > 1e-6:
        raise ValueError("Constrained expansion or work identity failed")
    return {"flux":individual, "expansion":trace}


def centered(positive, negative, scale):
    return {k: {f: (a-negative[k][f])/(2*scale) for f, a in v.items()}
            for k, v in positive.items()}


def report(folder):
    base, pairs = {}, {}
    for file in sorted(folder.glob("*.npz")):
        name = file.stem
        with np.load(file, allow_pickle=False) as data:
            if "Hraw" not in data or "z_inside" not in data:
                raise ValueError("Expected full-cadence native output: "+file.name)
            row = measure(data)
            lam = float(data["lam"])
        if lam == 0:
            base[name.removesuffix("_zero")] = row
        else:
            group = name.rsplit("_", 1)[0]
            pairs.setdefault(group, {}).setdefault(abs(lam), {})[np.sign(lam)] = row
    if not base:
        raise ValueError("No zero-coupling control")
    derived = {}
    for group, variations in pairs.items():
        derived[group] = {}
        if group.rsplit("_", 2)[0] not in base:
            raise ValueError("Missing matching zero-coupling control for " + group)
        for scale, signed in sorted(variations.items()):
            if set(signed) != {-1., 1.}:
                raise ValueError(f"Unpaired signed strengths for {group} at {scale}")
            derived[group][str(scale)] = centered(signed[1.], signed[-1.], scale)
    return {"proper_cut":CUT, "flux_units":"native kappa5^2 per normalized comoving volume",
            "zero_coupling":base, "derivative_per_y_squared":derived}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("runs", type=Path)
    ap.add_argument("--output", type=Path, default=HERE / "AUDIT.json")
    args = ap.parse_args()
    result = report(args.runs)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
