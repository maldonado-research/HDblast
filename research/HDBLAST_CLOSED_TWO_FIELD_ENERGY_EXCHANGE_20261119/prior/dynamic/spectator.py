"""Finite-band spectator propagation on the archived evolving closed slices.

The convention is the instantaneous Hamiltonian of the original daughter
field, in the archived unit-lapse geodesic proper time.  The output is a
finite-dimensional Bogoliubov diagnostic, not a radiation energy estimate.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("intrinsic_frequencies", HERE.parent/"intrinsic/frequencies.py")
intrinsic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(intrinsic)

ELL = 2
LOW_LABELS = (2, 4, 8)


class FieldHamiltonian:
    """Galerkin mass and stiffness in a fixed, physical radial basis.

    trial_basis is normalized with 2 sin(chi)^2 dchi/pi.  Multiplying its
    functions by sqrt(2/pi), and taking angular Y_lm of unit norm, gives the
    physical S3 normalization and exactly the matrices implemented below.
    """
    def __init__(self, geometry, solution, maximum=24, points=128, kinetic="matter"):
        self.geometry, self.solution = geometry, solution
        self.maximum, self.points, self.kinetic = maximum, points, kinetic
        self.labels = np.arange(ELL, maximum+1, 2)
        y = np.cos(np.pi*np.arange(1, points+1)/(points+1))
        self.y, self.sine2 = y, 1-y*y
        self.weights = 2*self.sine2/(points+1)
        self.interpolation = geometry.matrices(y)[0]
        v, d = intrinsic.trial_basis(ELL, maximum, y, self.weights)
        self.value, self.derivative = v[:, ::2], d[:, ::2]

    def matrices(self, t):
        values = self.interpolation @ self.solution.sol(t).reshape(6, self.geometry.n).T
        m, d, _, _, q, _ = values.T
        radial = m+self.sine2*d
        if self.kinetic == "matter":
            f = self.geometry.branch.fitted.spline(self.geometry.branch.top+q)[:, 1]
            z = 1/f
        elif self.kinetic == "geometric":
            z = np.ones_like(q)
        else:
            raise ValueError("Unknown kinetic convention")
        if not np.isfinite(z).all() or min(z) <= 0:
            raise ArithmeticError("Daughter kinetic weight left the positive fitted domain")
        a = self.geometry.background(t)[0]
        def gram(left, weight):
            return left.T @ (weight[:, None]*left)
        mass = a**3*gram(self.value, self.weights*z*np.exp(radial+2*m))
        stiffness = a*gram(self.derivative, self.weights*z*np.exp(2*m-radial))
        stiffness += a*ELL*(ELL+1)*gram(self.value, self.weights*z*np.exp(radial)/self.sine2)
        return mass, stiffness

    def modes(self, t):
        mass, stiffness = self.matrices(t)
        squared, vectors = eigh(stiffness, mass)
        if min(squared) <= 0:
            raise ArithmeticError("Finite-band Hamiltonian is not positive")
        # Reproducible real phases; no accidental degeneracy in this ell block.
        pivots = np.argmax(abs(vectors), axis=0)
        signs = np.sign(vectors[pivots, np.arange(len(pivots))])
        return mass, np.sqrt(squared), vectors*signs


def commutators(field, momentum):
    count = field.shape[0]
    xp = field @ momentum.conj().T-field.conj() @ momentum.T
    xx = field @ field.conj().T-field.conj() @ field.T
    pp = momentum @ momentum.conj().T-momentum.conj() @ momentum.T
    # Dimensionless versions avoid assigning significance to the field units.
    scales = [max(np.linalg.norm(field)**2, 1e-30), max(np.linalg.norm(momentum)**2, 1e-30)]
    return dict(xp_max_abs=float(np.max(abs(xp-1j*np.eye(count)))),
                xx_relative=float(np.max(abs(xx))/scales[0]),
                pp_relative=float(np.max(abs(pp))/scales[1]))


def initial_modes(hamiltonian):
    mass, omega, vectors = hamiltonian.modes(0.)
    field = vectors/np.sqrt(2*omega)[None, :]
    momentum = -1j*(mass @ vectors)*np.sqrt(omega/2)[None, :]
    return field, momentum


def project(hamiltonian, t, field, momentum):
    mass, omega, vectors = hamiltonian.modes(t)
    q = vectors.T @ mass @ field
    p = vectors.T @ momentum
    alpha = np.sqrt(omega/2)[:, None]*q+1j*p/np.sqrt(2*omega)[:, None]
    beta = np.sqrt(omega/2)[:, None]*q-1j*p/np.sqrt(2*omega)[:, None]
    return omega, alpha, beta


def propagate(hamiltonian, times, tolerance=2e-11, max_step=.006):
    field, momentum = initial_modes(hamiltonian)
    count = field.shape[0]
    initial = np.stack([field, momentum]).ravel()
    def rhs(t, state):
        field, momentum = state.reshape(2, count, count)
        mass, stiffness = hamiltonian.matrices(t)
        return np.stack([np.linalg.solve(mass, momentum), -stiffness @ field]).ravel()
    sol = solve_ivp(rhs, (0., max(times)), initial, method="DOP853",
                    rtol=tolerance, atol=2e-14, max_step=max_step, dense_output=True)
    if not sol.success:
        raise ArithmeticError(sol.message)
    rows = {}
    for t in times:
        u, p = sol.sol(t).reshape(2, count, count)
        omega, alpha, beta = project(hamiltonian, t, u, p)
        square = abs(beta)**2
        occupation = square.sum(axis=1)
        diag = np.diag(square)
        null = np.ones_like(square, dtype=bool)
        np.fill_diagonal(null, False)
        low = {str(n):dict(occupation=float(occupation[i]), diagonal=float(diag[i]),
                          from_other_initial_modes=float(occupation[i]-diag[i]))
               for i, n in enumerate(hamiltonian.labels) if n in LOW_LABELS}
        identity = np.eye(count)
        # These are column-mode identities, independent of particle conventions.
        bog = alpha.conj().T @ alpha-beta.conj().T @ beta-identity
        symmetric = alpha.T @ beta-beta.T @ alpha
        rows[f"{t:g}"] = dict(x=t, omega=omega.tolist(), low_modes=low,
                             beta_abs_squared=square.tolist(),
                             off_diagonal_beta_squared_sum=float(square[null].sum()),
                             off_diagonal_alpha_squared_sum=float((abs(alpha)**2)[null].sum()),
                             commutators=commutators(u, p),
                             bogoliubov_column_norm_error=float(np.max(abs(bog))),
                             bogoliubov_column_pair_error=float(np.max(abs(symmetric))))
    return rows, sol


def de_sitter_modes(h, labels, t, initial_field, initial_momentum):
    """Independent exact massless minimal solution in closed de Sitter.

    eta=atan(sinh(ht)), and chi_n proportional to
    (sin eta - i(n+1)cos eta) exp[-i(n+1)eta].
    The real F,G below enforce arbitrary initial chi and proper-time velocity.
    """
    eta = np.arctan(np.sinh(h*t))
    l = labels+1
    k2 = labels*(labels+2)
    c, s = np.cos(eta), np.sin(eta)
    phase = l*eta
    f = c*np.cos(phase)+(s/l)*np.sin(phase)
    g = (l*c*np.sin(phase)-s*np.cos(phase))/(h*k2)
    fdot = -h*(k2/l)*c*c*np.sin(phase)
    gdot = c*c*np.cos(phase)
    velocity = initial_momentum*h**3
    field = f[:, None]*initial_field+g[:, None]*velocity
    a = np.cosh(h*t)/h
    momentum = a**3*(fdot[:, None]*initial_field+gdot[:, None]*velocity)
    return field, momentum


def compare_rows(reference, other):
    return {str(n):dict(occupation_absolute_difference=float(other["low_modes"][str(n)]["occupation"]-
                                                   reference["low_modes"][str(n)]["occupation"]),
                        beta_row_max_abs_squared_difference=float(np.max(abs(
                            np.asarray(other["beta_abs_squared"])[i, :min(len(other["omega"]),len(reference["omega"]))]-
                            np.asarray(reference["beta_abs_squared"])[i, :min(len(other["omega"]),len(reference["omega"]))]))))
            for i, n in enumerate(range(ELL, 2*len(reference["omega"])+ELL, 2)) if n in LOW_LABELS}


def comparison_maximum(comparisons):
    return max(abs(item["occupation_absolute_difference"])
               for row in comparisons.values() for item in row.values())


def late_case(module, branch):
    """A resolved early-passage neighbourhood on the negative-minimum branch."""
    times = [0., .47]
    primary = module.Evolution(branch, 48, 192)
    finer = module.Evolution(branch, 80, 256)
    result = dict(x=.47, geometry_degree=48, refined_geometry_degree=80, runs={})
    for amplitude in (.001, -.001):
        record, solution = primary.evolve(amplitude, endpoint=.47, tolerance=2e-12, max_step=.003)
        fine_record, fine_solution = finer.evolve(amplitude, endpoint=.47, tolerance=2e-12, max_step=.003)
        entry = result["runs"][f"{amplitude:g}"] = dict(
            geometry_max_sampled_constraints=record["max_sampled_constraints"],
            refined_geometry_max_sampled_constraints=fine_record["max_sampled_constraints"],
            maximum_field=record["at_stop"]["maximum_field"], kinetic_controls={})
        for kinetic in ("matter", "geometric"):
            ham = FieldHamiltonian(primary, solution, 24, 160, kinetic)
            rows, ode = propagate(ham, times)
            item = entry["kinetic_controls"][kinetic] = dict(
                labels=ham.labels.tolist(), samples=rows, ode_evaluations=ode.nfev)
            item["band_refinements"] = []
            for cutoff in (12, 18, 30):
                other, _ = propagate(FieldHamiltonian(primary, solution, cutoff, 160, kinetic), times)
                item["band_refinements"].append(dict(maximum_n=cutoff,
                    comparisons={text:compare_rows(rows[text],other[text]) for text in rows}))
            fine, _ = propagate(FieldHamiltonian(finer, fine_solution, 24, 224, kinetic), times)
            item["geometry_refinement"] = {text:compare_rows(rows[text],fine[text]) for text in rows}
            time, _ = propagate(ham, times, tolerance=2e-12, max_step=.003)
            item["time_refinement"] = {text:compare_rows(rows[text],time[text]) for text in rows}
            quad, _ = propagate(FieldHamiltonian(primary, solution, 24, 224, kinetic), times)
            item["quadrature_refinement"] = {text:compare_rows(rows[text],quad[text]) for text in rows}
    return result


def verify(result):
    canonical, null, exact, refinement, bands, homogeneous_offdiagonal = [], [], [], [], [], []
    def inspect(run, homogeneous=False):
        for text, row in run["samples"].items():
            canonical.extend(row["commutators"].values())
            canonical.extend([row["bogoliubov_column_norm_error"],row["bogoliubov_column_pair_error"]])
            if text == "0":
                null.append(float(np.max(row["beta_abs_squared"])))
            if homogeneous:
                homogeneous_offdiagonal.extend([row["off_diagonal_beta_squared_sum"],
                                               row["off_diagonal_alpha_squared_sum"]])
        exact.extend(row["relative_exact_transfer_error"]
                     for row in run.get("closed_de_sitter_exact_control", []))
        for name in ("time_refinement", "quadrature_refinement", "geometry_refinement"):
            if name in run:
                refinement.append(comparison_maximum(run[name]))
        for band in run.get("band_refinements", []):
            if band["maximum_n"] >= 18:
                bands.append(comparison_maximum(band["comparisons"]))
    for entry in result["branches"].values():
        for amplitude, runs in entry["runs"].items():
            for run in runs.values():
                inspect(run, homogeneous=float(amplitude)==0.)
        if "early_passage_neighbourhood" in entry:
            for entry in entry["early_passage_neighbourhood"]["runs"].values():
                for run in entry["kinetic_controls"].values():
                    inspect(run)
    checks = dict(maximum_canonical_error=max(canonical),
                  maximum_null_time_occupation=max(null),
                  maximum_de_sitter_exact_transfer_relative_error=max(exact),
                  maximum_homogeneous_off_diagonal_squared_sum=max(homogeneous_offdiagonal),
                  maximum_time_quadrature_or_geometry_low_occupation_absolute_difference=max(refinement),
                  maximum_n18_or_higher_low_occupation_absolute_difference=max(bands))
    limits = dict(maximum_canonical_error=2e-9, maximum_null_time_occupation=1e-24,
                  maximum_de_sitter_exact_transfer_relative_error=2e-9,
                  maximum_homogeneous_off_diagonal_squared_sum=1e-24,
                  maximum_time_quadrature_or_geometry_low_occupation_absolute_difference=1e-9,
                  maximum_n18_or_higher_low_occupation_absolute_difference=1e-9)
    checks["required_limits"] = limits
    checks["passed"] = all(checks[name]<limit for name, limit in limits.items())
    if not checks["passed"]:
        raise ArithmeticError(f"Spectator control failed: {checks}")
    return checks


def calculate(gravity, extended=True):
    module = intrinsic.gravity_module(gravity)
    result = dict(convention="instantaneous original-field Hamiltonian in unit-lapse Einstein-frame geodesic proper time x",
                  kinetic_weight="B(Theta)=1/f(Theta), no derivative-free canonical field redefinition",
                  angular_ell=ELL, azimuthal_blocks="one m block; all 2ell+1 blocks identical",
                  radial_parity="even n-ell", primary_band_maximum=24,
                  occupation="sum over initial positive-frequency basis columns, within the stated finite band",
                  interpretation="finite-time, finite-band Bogoliubov diagnostic; no UV sum, thermal energy or backreaction",
                  branches={})
    for key, branch in module.constraint.load_branches(gravity/"previous").items():
        geometry = module.Evolution(branch, 32, 128)
        entry = result["branches"][key] = dict(runs={})
        for amplitude in (0., .001, -.001):
            times = [0., .12, .20] if extended else [0., .12]
            _, solution = geometry.evolve(amplitude, endpoint=max(times))
            runs = entry["runs"][f"{amplitude:g}"] = {}
            for kinetic in ("matter", "geometric"):
                hamiltonian = FieldHamiltonian(geometry, solution, 24, 128, kinetic)
                rows, ode = propagate(hamiltonian, times)
                runs[kinetic] = dict(labels=hamiltonian.labels.tolist(), samples=rows,
                                     ode_evaluations=ode.nfev)
                if amplitude == 0:
                    u0, p0 = initial_modes(hamiltonian)
                    exact = []
                    for t in times:
                        # Z is a constant in this case, so its sole effect is
                        # P=Z a^3 dot X. Correct the analytic velocity mapping.
                        z = 1/branch.fitted.spline(branch.top)[1] if kinetic == "matter" else 1.
                        u, p = de_sitter_modes(branch.h, hamiltonian.labels, t, u0, p0/z)
                        p *= z
                        un, pn = ode.sol(t).reshape(2, len(hamiltonian.labels), -1)
                        error = max(float(np.max(abs(un-u))/np.max(abs(u))),
                                    float(np.max(abs(pn-p))/np.max(abs(p))))
                        exact.append(dict(x=t, relative_exact_transfer_error=error))
                    runs[kinetic]["closed_de_sitter_exact_control"] = exact
                else:
                    # Finite-band convergence is evaluated for the low outputs,
                    # retaining each truncation's entire initial vacuum basis.
                    band = []
                    for cutoff in (12, 18):
                        other = FieldHamiltonian(geometry, solution, cutoff, 128, kinetic)
                        coarse, _ = propagate(other, times)
                        band.append(dict(maximum_n=cutoff,
                            comparisons={text:compare_rows(rows[text],coarse[text]) for text in rows}))
                    runs[kinetic]["band_refinements"] = band
                    temporal, _ = propagate(hamiltonian, times, tolerance=2e-12, max_step=.003)
                    runs[kinetic]["time_refinement"] = {text:compare_rows(rows[text],temporal[text]) for text in rows}
                    spatial_ham = FieldHamiltonian(geometry, solution, 24, 192, kinetic)
                    spatial, _ = propagate(spatial_ham, times)
                    runs[kinetic]["quadrature_refinement"] = {text:compare_rows(rows[text],spatial[text]) for text in rows}
            print(key, f"{amplitude:g}", "spectator complete", flush=True)
        # One geometry refinement for both signs and both kinetic choices.
        finer = module.Evolution(branch, 48, 192)
        for amplitude in (.001, -.001):
            _, solution = finer.evolve(amplitude, endpoint=max(times), tolerance=2e-12, max_step=.003)
            for kinetic in ("matter", "geometric"):
                rows, _ = propagate(FieldHamiltonian(finer, solution, 24, 192, kinetic), times)
                primary = entry["runs"][f"{amplitude:g}"][kinetic]
                primary["geometry_refinement"] = {text:compare_rows(primary["samples"][text],rows[text]) for text in rows}
        if extended and branch.fitted.pm < 0:
            entry["early_passage_neighbourhood"] = late_case(module, branch)
    result["checks"] = verify(result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gravity", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--short", action="store_true", help="Only the x=0.12 primary endpoint")
    args = parser.parse_args()
    data = calculate(args.gravity.resolve(), extended=not args.short)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
