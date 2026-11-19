"""Independent real-state extremal signed coupling work across first passage.

Only the archived background and radial trial basis are imported. Daughter
transfer, split work forms, initial energy whitening and eigenvalue tests are
implemented here. No vacuum subtraction or particle convention is needed.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

HERE = Path(__file__).resolve().parent


def load(path):
    spec = importlib.util.spec_from_file_location("phase_check_"+path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def positive_power(matrix, power):
    values, vectors = eigh(matrix)
    if values.min() <= 0:
        raise ArithmeticError("Initial energy is not positive")
    return (vectors*values**power) @ vectors.T


class WorkHamiltonian:
    """Original field q and canonical p on one real ell=2 angular block."""
    def __init__(self, passage, geometry, solution, maximum, points, geometric=False):
        self.geometry, self.solution, self.geometric = geometry, solution, geometric
        self.labels = np.arange(2, maximum+1, 2)
        self.y = np.cos(np.pi*np.arange(1, points+1)/(points+1))
        self.w = 1-self.y**2
        self.weights = 2*self.w/(points+1)
        self.interpolation = geometry.matrices(self.y)[0]
        value, derivative = passage.dynamic.intrinsic.trial_basis(
            2, maximum, self.y, self.weights)
        self.value, self.derivative = value[:, ::2], derivative[:, ::2]
        self.f_top = float(geometry.branch.fitted.spline(geometry.branch.top)[1])

    def matrices(self, t):
        g = self.geometry
        m, d, b, e, q, pi = (self.interpolation @
                             self.solution.sol(t).reshape(6, g.n).T).T
        a, H, _ = g.background(t)
        if self.geometric:
            z, rate = np.ones_like(q), np.zeros_like(q)
        else:
            spline = g.branch.fitted.spline
            f = spline(g.branch.top+q)[:, 1]
            z = self.f_top/f
            rate = -spline(g.branch.top+q, 1)[:, 1]*pi/f
        if not np.isfinite(z).all() or z.min() <= 0:
            raise ArithmeticError("Positive fitted kinetic domain left")
        volume = a**3*self.weights*z*np.exp(3*m+self.w*d)
        radial = a*self.weights*z*np.exp(m-self.w*d)
        angular = a*6*self.weights*z*np.exp(m+self.w*d)/self.w

        def gram(left, weight):
            return left.T @ (weight[:, None]*left)

        mass = gram(self.value, volume)
        stiff = gram(self.derivative, radial)+gram(self.value, angular)
        mz = gram(self.value, volume*rate)
        kz = gram(self.derivative, radial*rate)+gram(self.value, angular*rate)
        mg = gram(self.value, volume*(3*H+3*b+self.w*e))
        kg = gram(self.derivative, radial*(H+b-self.w*e))
        kg += gram(self.value, angular*(H+b+self.w*e))
        return mass, stiff, mz, kz, mg, kg


def energy_form(ham, t, transfer):
    mass, stiff, *_ = ham.matrices(t)
    count = len(ham.labels)
    q, p = transfer[:count], transfer[count:]
    return q.T @ stiff @ q+p.T @ np.linalg.solve(mass, p)


def summarize_form(work, before, after, initial_map, labels):
    work = (work+work.T)/2
    values, vectors = eigh(work)
    bracket_values, bracket_vectors = eigh(work, (before+before.T)/2)

    def state_record(phase, name):
        state = initial_map @ phase
        # One shared real eigenvector phase; stabilize replay across LAPACKs.
        pivot = int(np.argmax(abs(state)))
        if state[pivot] < 0:
            phase, state = -phase, -state
        count = len(labels)
        return dict(extremum=name, initial_energy=float(phase @ phase/2),
            coupling_work=float(phase @ work @ phase/2),
            bracket_start_energy=float(phase @ before @ phase/2),
            bracket_end_energy=float(phase @ after @ phase/2),
            initial_field_coefficients=state[:count].tolist(),
            initial_canonical_momentum_coefficients=state[count:].tolist())

    # First pair has E(bounce)=1; second pair has E(bracket start)=1.
    states = [state_record(np.sqrt(2)*vectors[:, position], name)
              for position, name in ((0, "minimum"), (-1, "maximum"))]
    bracket_states = [state_record(np.sqrt(2)*bracket_vectors[:, position], name)
                     for position, name in ((0, "minimum"), (-1, "maximum"))]
    return dict(bounce_energy_ratio_minimum=float(values[0]),
        bounce_energy_ratio_maximum=float(values[-1]),
        bracket_start_energy_ratio_minimum=float(bracket_values[0]),
        bracket_start_energy_ratio_maximum=float(bracket_values[-1]),
        eigenvalues_over_bounce_energy=values.tolist(), extremizing_states=states,
        bracket_start_extremizing_states=bracket_states)


def run(ham, times, initial_maximum=None, tolerance=2e-11, max_step=.006,
        witness=None):
    """Propagate the whole real transfer; accumulate work only within bracket."""
    mass, stiff, *_ = ham.matrices(0.)
    count = len(ham.labels)
    low = count if initial_maximum is None else int(np.sum(ham.labels <= initial_maximum))
    embed = np.eye(count)[:, :low]
    # Specify fixed initial support in field and velocity, rather than momentum.
    # Thus a support n<=8 means the same physical initial data at larger cutoffs.
    initial_map = np.zeros((2*count, 2*low))
    initial_map[:count, :low] = embed @ positive_power(stiff[:low, :low], -.5)
    initial_map[count:, low:] = mass @ embed @ positive_power(mass[:low, :low], -.5)
    dimension = 2*low

    def derivative(t, transfer, work=False):
        mass, stiff, mz, kz, mg, kg = ham.matrices(t)
        q, p = transfer[:count], transfer[count:]
        velocity = np.linalg.solve(mass, p)
        evolved = np.vstack((velocity, -stiff @ q))
        if not work:
            return evolved
        coupling = q.T @ kz @ q-velocity.T @ mz @ velocity
        geometric = q.T @ kg @ q-velocity.T @ mg @ velocity
        return evolved, coupling, geometric

    before_t, after_t = times
    first = solve_ivp(lambda t, y: derivative(t, y.reshape(2*count, dimension)).ravel(),
        (0., before_t), initial_map.ravel(), method="DOP853", rtol=tolerance,
        atol=2e-14, max_step=max_step)
    if not first.success:
        raise ArithmeticError(first.message)
    before_transfer = first.y[:, -1].reshape(2*count, dimension)
    before = energy_form(ham, before_t, before_transfer)
    size = 2*count*dimension
    initial = np.r_[before_transfer.ravel(), np.zeros(2*dimension**2)]

    def rhs(t, state):
        evolved, coupling, geometric = derivative(
            t, state[:size].reshape(2*count, dimension), True)
        return np.r_[evolved.ravel(), coupling.ravel(), geometric.ravel()]

    second = solve_ivp(rhs, times, initial, method="DOP853", rtol=tolerance,
                       atol=2e-14, max_step=max_step/4)
    if not second.success:
        raise ArithmeticError(second.message)
    state = second.y[:, -1]
    after_transfer = state[:size].reshape(2*count, dimension)
    wz, wg = state[size:].reshape(2, dimension, dimension)
    after = energy_form(ham, after_t, after_transfer)
    identity_error = float(np.max(abs(wz+wg-(after-before))))
    initial_energy_error = float(np.max(abs(energy_form(ham, 0., initial_map)-np.eye(dimension))))
    result = summarize_form(wz, before, after, initial_map, ham.labels)
    result.update(propagated_labels=ham.labels.tolist(),
        initial_field_and_velocity_maximum_n=int(ham.labels[low-1]),
        ode_evaluations=first.nfev+second.nfev,
        ledger_matrix_max_absolute_error=identity_error,
        initial_energy_whitening_max_absolute_error=initial_energy_error,
        geometry_work_matrix_eigenvalues=eigh((wg+wg.T)/2, eigvals_only=True).tolist(),
        maximum_coupling_matrix_asymmetry=float(np.max(abs(wz-wz.T))))
    if witness is not None:
        phase = np.r_[positive_power(stiff[:low, :low], .5) @ witness["field"],
                      positive_power(mass[:low, :low], .5) @ witness["velocity"]]
        e0 = float(phase @ phase/2)
        eb = float(phase @ before @ phase/2)
        work = float(phase @ wz @ phase/2)
        result["unchanged_n24_bracket_maximizing_witness"] = dict(
            initial_bounce_energy=e0, bracket_start_energy=eb,
            bracket_end_energy=float(phase @ after @ phase/2), coupling_work=work,
            work_over_bounce_energy=work/e0, work_over_bracket_start_energy=work/eb)
    return result


def endpoint_difference(reference, other):
    names = ("bounce_energy_ratio_minimum", "bounce_energy_ratio_maximum",
             "bracket_start_energy_ratio_minimum", "bracket_start_energy_ratio_maximum")
    return {name:other[name]-reference[name] for name in names}


def calculate(prior):
    passage = load(prior/"passage/passage.py")
    gravity = prior/"gravity"
    module = passage.dynamic.intrinsic.gravity_module(gravity)
    key, branch = next((key, branch) for key, branch in
        module.constraint.load_branches(gravity/"previous").items() if branch.fitted.pm < 0)
    result = dict(branch=key, angular_ell=2, one_real_angular_block=True,
        geometry_normalization="canonical scalar action, V=3u; initial energy reference 3h^2",
        kinetic_weight="Z=f(Theta_top)/f(Theta), positive constant normalization of archived inverse f",
        work_sign="positive is net coupling energy into daughter; geometry work excluded",
        initial_energy="positive original-field Hamiltonian on t=0 bounce slice; bracket-start ratios also supplied",
        bracket_half_width=passage.HALF_WIDTH,
        statement="finite radial bands of arbitrary real linear states, not a UV sum or radiation fraction",
        runs={})
    for amplitude in (.001, -.001):
        geometry, solution, event = passage.geometry_through_passage(
            module, branch, amplitude, 64, 160)
        times = (event["sample_times"]["before"], event["sample_times"]["after"])
        entry = dict(crossing_x=event["crossing_x"], bracket=list(times),
            maximum_sampled_constraints=event["maximum_sampled_constraints"],
            full_initial_bands={}, fixed_low_support={})
        for cutoff in (8, 12, 18, 24, 36):
            ham = WorkHamiltonian(passage, geometry, solution, cutoff, 192)
            if cutoff <= 24:
                entry["full_initial_bands"][str(cutoff)] = run(ham, times)
            if cutoff == 8:
                entry["fixed_low_support"][str(cutoff)] = entry["full_initial_bands"][str(cutoff)]
            else:
                entry["fixed_low_support"][str(cutoff)] = run(ham, times, initial_maximum=8)
        reference = entry["full_initial_bands"]["18"]
        ham = WorkHamiltonian(passage, geometry, solution, 18, 192)
        entry["time_refinement"] = endpoint_difference(reference,
            run(ham, times, tolerance=2e-12, max_step=.003))
        entry["quadrature_refinement"] = endpoint_difference(reference,
            run(WorkHamiltonian(passage, geometry, solution, 18, 272), times))
        finer, fine_solution, fine_event = passage.geometry_through_passage(
            module, branch, amplitude, 96, 256, tolerance=2e-13, max_step=.0015)
        fine_times = (fine_event["sample_times"]["before"], fine_event["sample_times"]["after"])
        entry["geometry_refinement"] = endpoint_difference(reference,
            run(WorkHamiltonian(passage, finer, fine_solution, 18, 272), fine_times))
        entry["event_time_refinement_absolute_difference"] = abs(fine_event["crossing_x"]-event["crossing_x"])
        entry["constant_Z_control"] = run(
            WorkHamiltonian(passage, geometry, solution, 18, 192, geometric=True), times)
        entry["fixed_support_n24_to_n36_refinement"] = endpoint_difference(
            entry["fixed_low_support"]["24"], entry["fixed_low_support"]["36"])
        # Hold the near-envelope n<=24 initial data fixed at larger dynamical cutoffs.
        low_ham = WorkHamiltonian(passage, geometry, solution, 24, 192)
        mass24 = low_ham.matrices(0.)[0]
        maximizing = entry["full_initial_bands"]["24"]["bracket_start_extremizing_states"][-1]
        witness = dict(field=np.asarray(maximizing["initial_field_coefficients"]),
            velocity=np.linalg.solve(mass24,
                np.asarray(maximizing["initial_canonical_momentum_coefficients"])))
        support24 = entry["fixed_n24_support_high_band_check"] = dict(
            initial_field_maximum_n=24, initial_velocity_maximum_n=24,
            unchanged_witness_field_coefficients=witness["field"].tolist(),
            unchanged_witness_velocity_coefficients=witness["velocity"].tolist(),
            original_n24_witness_work_over_bracket_start_energy=float(
                maximizing["coupling_work"]/maximizing["bracket_start_energy"]),
            propagated_bands={})
        for cutoff in (36, 48):
            support24["propagated_bands"][str(cutoff)] = run(
                WorkHamiltonian(passage, geometry, solution, cutoff, 192), times,
                initial_maximum=24, witness=witness)
        support24["n24_to_n36_extremum_difference"] = endpoint_difference(
            entry["full_initial_bands"]["24"], support24["propagated_bands"]["36"])
        support24["n36_to_n48_extremum_difference"] = endpoint_difference(
            support24["propagated_bands"]["36"], support24["propagated_bands"]["48"])
        witness36, witness48 = [support24["propagated_bands"][str(cutoff)][
            "unchanged_n24_bracket_maximizing_witness"] for cutoff in (36, 48)]
        support24["n36_to_n48_unchanged_witness_bracket_ratio_difference"] = (
            witness48["work_over_bracket_start_energy"]-
            witness36["work_over_bracket_start_energy"])
        support24["n24_to_n48_unchanged_witness_bracket_ratio_difference"] = (
            witness48["work_over_bracket_start_energy"]-
            support24["original_n24_witness_work_over_bracket_start_energy"])
        result["runs"][f"{amplitude:g}"] = entry
        print(f"seed {amplitude:g}: full and fixed-support work bounds complete", flush=True)
    results = [row for entry in result["runs"].values()
               for group in ("full_initial_bands", "fixed_low_support")
               for row in entry[group].values()]
    results.extend(row for entry in result["runs"].values() for row in
        entry["fixed_n24_support_high_band_check"]["propagated_bands"].values())
    refinements = [abs(v) for entry in result["runs"].values()
        for name in ("time_refinement", "quadrature_refinement", "geometry_refinement")
        for v in entry[name].values()]
    controls = dict(maximum_ledger_matrix_absolute_error=max(
        row["ledger_matrix_max_absolute_error"] for row in results),
        maximum_initial_energy_whitening_error=max(
        row["initial_energy_whitening_max_absolute_error"] for row in results),
        maximum_refinement_extremum_absolute_difference=max(refinements),
        maximum_event_time_refinement_absolute_difference=max(
        entry["event_time_refinement_absolute_difference"] for entry in result["runs"].values()),
        maximum_fixed_support_n24_to_n36_extremum_absolute_difference=max(
        abs(value) for entry in result["runs"].values()
        for value in entry["fixed_support_n24_to_n36_refinement"].values()),
        maximum_n24_support_n36_to_n48_extremum_absolute_difference=max(
        abs(value) for entry in result["runs"].values() for value in
        entry["fixed_n24_support_high_band_check"]["n36_to_n48_extremum_difference"].values()),
        maximum_n24_support_n36_to_n48_unchanged_witness_bracket_ratio_absolute_difference=max(
        abs(entry["fixed_n24_support_high_band_check"][
            "n36_to_n48_unchanged_witness_bracket_ratio_difference"])
        for entry in result["runs"].values()),
        maximum_constant_Z_coupling_absolute_extremum=max(
        abs(entry["constant_Z_control"][name]) for entry in result["runs"].values()
        for name in ("bounce_energy_ratio_minimum", "bounce_energy_ratio_maximum")),
        both_signs_realized_in_every_tested_band=all(
        row["bounce_energy_ratio_minimum"] < 0 < row["bounce_energy_ratio_maximum"] for row in results))
    controls["passed"] = (controls["maximum_ledger_matrix_absolute_error"] < 2e-9
        and controls["maximum_initial_energy_whitening_error"] < 2e-11
        and controls["maximum_refinement_extremum_absolute_difference"] < 2e-8
        and controls["maximum_event_time_refinement_absolute_difference"] < 2e-9
        and controls["maximum_fixed_support_n24_to_n36_extremum_absolute_difference"] < 2e-9
        and controls["maximum_n24_support_n36_to_n48_extremum_absolute_difference"] < 2e-8
        and controls["maximum_n24_support_n36_to_n48_unchanged_witness_bracket_ratio_absolute_difference"] < 2e-8
        and controls["maximum_constant_Z_coupling_absolute_extremum"] == 0.
        and controls["both_signs_realized_in_every_tested_band"])
    result["checks"] = controls
    if not controls["passed"]:
        raise ArithmeticError(f"Independent phase work checks failed: {controls}")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, default=HERE.parent/"prior")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = calculate(args.prior.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
