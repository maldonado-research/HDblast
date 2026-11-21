"""Two canonical scalar fluctuation fields on the HDblast #192 short-drive background.

Run with Python 3, NumPy and SciPy. All momenta and frequencies use m as unit.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp

MASS = 1.9
AMPLITUDE = .008
WIDTH = .24
INCOMING = .12287824489895026
REFERENCE_PARENT = .00027167585541013783
EDGES = np.array([1 / 8, 1 / 4, 1 / 2, 1., 2., 4., 8., 16., 32., 64.])
# ENVELOPE.json/inner_source and SWEEP.json/runs at HDblast dc99ff5.
ACTIONS = {"4e10": 3.608589582002018e48, "4e12": 1.0415592372880226e54}
SWEEP_CAPS = {"4e10": 1.2659228153859582e46, "4e12": 1.2659138001404031e54}
PRIMARY_STEPS = 1024
PRIMARY_PANELS = 128


def step(x):
    """Smooth sigmoid with its first two derivatives. Use the small exponential at both ends."""
    if not 0 < x < 1:
        return float(x >= 1), 0., 0.
    y = 1 - x
    z = 1 / y - 1 / x
    if abs(z) > 740:
        return float(z > 0), 0., 0.
    e = math.exp(-abs(z))
    b = (1 if z >= 0 else e) / (1 + e)
    d = e / (1 + e)**2
    dz = 1 / x**2 + 1 / y**2
    ddz = 2 / y**3 - 2 / x**3
    return b, d * dz, d * (ddz + (1 - 2 * b) * dz**2)


def bump(u):
    scale = 2 / WIDTH
    x = u * scale
    up = step(x)
    down = step(x - 1)
    return tuple((a - b) * scale**n for n, (a, b) in enumerate(zip(up, down)))


def background():
    def rhs(s, state):
        u, c, velocity, momentum = state
        b, db, _ = bump(u)
        force = MASS**2 * (c - AMPLITUDE * b)
        return velocity, momentum, force * AMPLITUDE * db, -force

    def exit_drive(s, state):
        return state[0] - WIDTH
    exit_drive.terminal = True
    exit_drive.direction = 1
    sol = solve_ivp(rhs, (0, 3), (0, 0, INCOMING, 0), method="DOP853",
                    rtol=2e-13, atol=2e-15, dense_output=True,
                    events=exit_drive, max_step=.02)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError("Background did not exit the drive")
    return sol


def hessian_perturbation(state):
    u, c = state[:2]
    b, db, ddb = bump(u)
    slope = AMPLITUDE * db
    return MASS**2 * np.array([[slope**2 - (c - AMPLITUDE*b)*AMPLITUDE*ddb, -slope],
                               [-slope, 0.]])


def evolve_modes(sol, k, steps):
    """Evolve interaction-picture alpha and beta by RK4, with frequencies exact in the phases."""
    omega = np.column_stack((k, np.hypot(k, MASS)))
    norm = 1 / np.sqrt(2 * omega)
    alpha = np.broadcast_to(np.eye(2), (len(k), 2, 2)).astype(complex).copy()
    beta = np.zeros_like(alpha)
    h = sol.t[-1] / steps
    times = np.linspace(0., sol.t[-1], 2 * steps + 1)
    potentials = [hessian_perturbation(state) for state in sol.sol(times).T]

    def rhs(index, a, b):
        neg = (np.exp(-1j * omega * times[index]) * norm)[:, :, None]
        pos = neg.conj()
        field = neg * a + pos * b
        force = potentials[index] @ field
        return -1j * pos * force, 1j * neg * force

    for index in range(0, 2 * steps, 2):
        a1, b1 = rhs(index, alpha, beta)
        a2, b2 = rhs(index + 1, alpha + h*a1/2, beta + h*b1/2)
        a3, b3 = rhs(index + 1, alpha + h*a2/2, beta + h*b2/2)
        a4, b4 = rhs(index + 2, alpha + h*a3, beta + h*b3)
        alpha += h/6 * (a1 + 2*a2 + 2*a3 + a4)
        beta += h/6 * (b1 + 2*b2 + 2*b3 + b4)
    return alpha, beta


def direct(sol, q):
    """Independent DOP853 evolution of the original position and velocity matrices."""
    omega = np.array([q, np.hypot(q, MASS)])
    init = np.diag(1 / np.sqrt(2 * omega))
    def rhs(t, vector):
        field, velocity = vector.reshape(2, 2, 2)
        freq2 = np.diag(omega**2) + hessian_perturbation(sol.sol(t))
        return np.stack((velocity, -freq2 @ field)).ravel()

    state = np.stack((init, -1j*omega[:, None]*init)).ravel()
    evolved = solve_ivp(rhs, (0., sol.t[-1]), state, method="DOP853",
                        rtol=2e-12, atol=2e-14, max_step=.013)
    if not evolved.success:
        raise RuntimeError(evolved.message)
    field, velocity = evolved.y[:, -1].reshape(2, 2, 2)
    return ((np.exp(-1j*omega*sol.t[-1]) / np.sqrt(2*omega))[:, None]
            * (omega[:, None]*field - 1j*velocity))


def calculation(sol, steps, panels):
    segments = [np.linspace(left, right, panels+1) for left, right in zip(EDGES[:-1], EDGES[1:])]
    k = np.unique(np.concatenate(segments))
    alpha, beta = evolve_modes(sol, k, steps)
    at, bt = alpha.transpose(0, 2, 1), beta.transpose(0, 2, 1)
    wronskian = np.max(abs(at.conj() @ alpha - bt.conj() @ beta - np.eye(2)))
    symplectic = np.max(abs(at @ beta - bt @ alpha))
    occ = np.sum(abs(beta)**2, axis=2)
    omega = np.column_stack((k, np.hypot(k, MASS)))
    integrand = (k**2 / (2*math.pi**2))[:, None] * omega * occ
    integrals = [simpson(integrand[i*panels:(i+1)*panels+1], x=s, axis=0)
                 for i, s in enumerate(segments)]
    return k, beta, np.array(integrals), float(wronskian), float(symplectic)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    sol = background()
    _, c, _, momentum = sol.y[:, -1]
    parent = .5 * (momentum**2 + (MASS*c)**2)
    states = sol.sol(np.linspace(0., sol.t[-1], 2049))
    total = .5 * (states[2]**2 + states[3]**2 +
                   (MASS * (states[1] - AMPLITUDE*np.array([bump(u)[0] for u in states[0]])))**2)
    k, beta, energies, wronskian, symplectic = calculation(sol, PRIMARY_STEPS, PRIMARY_PANELS)
    _, _, rough_time, _, _ = calculation(sol, PRIMARY_STEPS // 2, PRIMARY_PANELS)
    _, _, rough_k, _, _ = calculation(sol, PRIMARY_STEPS, PRIMARY_PANELS // 2)
    energy = np.sum(energies)
    checks = {}
    samples = {}
    for q in EDGES:
        matrix = beta[np.searchsorted(k, q)]
        direct_matrix = direct(sol, q)
        checks[str(q)] = float(np.max(abs(matrix-direct_matrix)))
        samples[str(q)] = (abs(matrix)**2).tolist()

    # Canonical normalization can be checked starting either in wall units
    # (energy density in sigma = S0 / L^4) or directly in fast quantum units.
    def energy_ratios(actions):
        ratios = {}
        for mass, action in actions.items():
            mu = float(mass)
            wall_answer = (mu**4 * energy) / (action * parent)
            canonical_background_energy_in_m4 = (math.sqrt(action) / mu**2)**2 * parent
            fast_answer = energy / canonical_background_energy_in_m4
            if not math.isclose(wall_answer, fast_answer, rel_tol=1e-14):
                raise ArithmeticError("Canonical and wall units disagree")
            ratios[mass] = wall_answer
        return ratios

    result = {
        "time_steps": PRIMARY_STEPS, "momentum_panels_per_interval": PRIMARY_PANELS,
        "incoming_velocity": INCOMING, "duration_in_fast_time": sol.t[-1],
        "classical_outgoing_parent_energy": parent,
        "classical_energy_max_relative_drift": float(np.max(abs(total/total[0]-1))),
        "parent_relative_to_pinned_minus_one": parent / REFERENCE_PARENT - 1,
        "momentum_edges_in_daughter_mass_units": EDGES.tolist(),
        "interval_energy_per_m4_by_outgoing_species_driver_parent": energies.tolist(),
        "interval_energy_sum_per_m4": energy,
        "cutoff_cumulative_sums_per_m4": np.cumsum(energies.sum(axis=1)).tolist(),
        "sample_abs_beta_squared_rows_out_columns_in": samples,
        "absolute_energy_sum_change_halving_time_steps": abs(energy - rough_time.sum()),
        "absolute_energy_sum_change_halving_momentum_panels": abs(energy - rough_k.sum()),
        "max_interval_energy_absolute_change_halving_time_steps": float(np.max(abs(energies-rough_time))),
        "max_interval_energy_absolute_change_halving_momentum_panels": float(np.max(abs(energies-rough_k))),
        "max_wronskian_error": wronskian, "max_bosonic_symmetry_error": symplectic,
        "direct_beta_max_absolute_difference_per_sample": checks,
        "band_energy_to_classical_parent": energy_ratios(ACTIONS),
        "band_energy_to_classical_parent_at_half_step_aggregate_sweep_ceiling": energy_ratios(SWEEP_CAPS),
        "pinned_half_step_aggregate_sweep_action_ceilings": SWEEP_CAPS,
        "actions": ACTIONS,
        "action_over_mu4_for_equal_band_and_classical_parent_energy": energy / parent,
    }
    if (abs(parent / REFERENCE_PARENT - 1) > 1e-8 or wronskian > 1e-7 or symplectic > 1e-7
            or max(checks.values()) > 1e-7 or max(abs(energy - rough_time.sum()),
                                                    abs(energy - rough_k.sum())) > 1e-5*energy):
        raise RuntimeError("Convergence or normalization check failed: " + repr(result))
    receipt = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(receipt)
    else:
        print(receipt, end="")


if __name__ == "__main__":
    main()
