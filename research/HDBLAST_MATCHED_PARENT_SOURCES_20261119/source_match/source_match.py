"""Finite-grid primal and dual bounds with initial pressure and scalar-force matching."""
import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import block_diag, csr_matrix, hstack, vstack


PARENT = Path("research/HDBLAST_MOVING_PARENT_20261119")
PINS = {"phase_feasibility/feasibility.py": "4b4196bb9316fa12fa6407da7c27df2b4e036486f627c3c74db93e5f7fe43bd7",
        "coherence/coherence.py": "9e94c224d22774c1073d0bac589bbcb50f1106c2f0d77e801c5fecab8f67a29a"}
LP_OPTIONS = {"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10}


def load_parent(repo):
    for filename, digest in PINS.items():
        path = repo / PARENT / filename
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected historical response source: {filename}")
    path = repo / PARENT / "phase_feasibility/feasibility.py"
    spec = importlib.util.spec_from_file_location("phase_response", path)
    parent = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parent)
    return parent


def initial_states(prior, repo):
    initial = {}
    with zipfile.ZipFile(io.BytesIO(prior.checked(repo / prior.WALL, prior.WALL_SHA))) as archive:
        for grid in ("fine", "coarse"):
            reference = json.loads(prior.checked(repo / prior.REF.format(grid), prior.REF_SHA[grid]))
            scale = reference["matter"]["events"][0]["A_star"]
            for arm in ("quiet", "plus", "minus", "far"):
                geom, start = prior.wall(archive, grid, arm, scale)
                a, _, mass = geom(start)
                state = dict(clock=float(start), a=float(a), mass=float(mass))
                if arm != "quiet" and initial[grid] != state:
                    raise ValueError(f"Different initial state on {arm}")
                initial[grid] = state
    return initial


def source_constraints(q, weights, initial, both_meshes):
    energies = {grid: np.hypot(q, row["a"] * row["mass"]) for grid, row in initial.items()}
    u, v = energies["fine"], energies["coarse"]
    physical = [weights/u, weights*q*q/u]
    if both_meshes:
        physical.extend((weights/v, weights*q*q/v))
        s = u + v
        # Exactly equivalent span for unequal initial (am)^2. Direct subtraction
        # of the nearly equal physical functions loses the highest moment.
        if initial["fine"]["a"]*initial["fine"]["mass"] == initial["coarse"]["a"]*initial["coarse"]["mass"]:
            raise ValueError("Use only the two independent constraints for identical meshes")
        work = weights * np.array([s, 1/s, 1/(u*v*s), 1/(u*v*s**3)])
    else:
        work = np.array(physical)
    scaled = work / np.max(abs(work), axis=1, keepdims=True)
    condition = float(np.linalg.cond(scaled))
    orthogonal = np.linalg.qr(scaled.T, mode="reduced")[0].T
    physical = np.asarray(physical)
    physical /= abs(physical).sum(axis=1, keepdims=True)
    return orthogonal, physical, condition


def bounds(matrix, moments=None, vertices=256):
    rows, nodes = matrix.shape
    imaginary = moments is None
    if imaginary:
        # z=i*b, so each signed margin is -Im(V)*b.
        inequalities = np.column_stack((matrix.imag, np.ones(rows)))
        result = linprog(np.r_[np.zeros(nodes), -1.], A_ub=inequalities, b_ub=np.zeros(rows),
                         bounds=[(-1., 1.)]*nodes + [(None, None)], method="highs", options=LP_OPTIONS)
    else:
        g, physical, condition = moments
        coefficients = np.stack((matrix.real, -matrix.imag), axis=-1).reshape(rows, 2*nodes)
        reversal = hstack((csr_matrix(-coefficients), csr_matrix(np.ones((rows, 1)))))
        angles = (2*np.arange(vertices) + 1)*np.pi/vertices
        sides = np.column_stack((np.cos(angles), np.sin(angles)))
        disks = hstack((block_diag([sides]*nodes), csr_matrix((nodes*vertices, 1))))
        equalities = np.zeros((len(g), 2*nodes + 1))
        equalities[:, :-1:2] = g
        result = linprog(np.r_[np.zeros(2*nodes), -1.], A_ub=vstack((reversal, disks), format="csr"),
                         b_ub=np.r_[np.zeros(rows), np.full(nodes*vertices, np.cos(np.pi/vertices))],
                         A_eq=equalities, b_eq=np.zeros(len(g)), bounds=(None, None),
                         method="highs", options=LP_OPTIONS)
    if not result.success:
        raise ValueError(result.message)
    p = np.maximum(-result.ineqlin.marginals[:rows], 0)
    p /= p.sum()
    summed = p @ matrix
    if imaginary:
        z = 1j * np.clip(result.x[:-1], -1., 1.)
        upper = float(abs(summed.imag).sum())
        multipliers = []
        residual = 0.
        basis_residual = 0.
        condition = None
    else:
        z = result.x[:-1:2] + 1j*result.x[1:-1:2]
        z -= g.T @ (g @ z.real)
        z /= max(1., abs(z).max())
        multipliers = result.eqlin.marginals
        upper = float(abs(summed + multipliers @ g).sum())
        residual = float(abs(physical @ z.real).max())
        basis_residual = float(abs(g @ z.real).max())
        multipliers = multipliers.tolist()
    margins = (matrix @ z).real
    lower = float(min(margins))
    if max(residual, basis_residual) > 1e-9 or max(abs(z)) > 1 + 1e-12 or lower > upper + 1e-10:
        raise ValueError("Primal or dual certificate failed")
    return dict(lower=lower, upper=upper, quantum_lower=float(np.sqrt(2)*lower),
                quantum_upper=float(np.sqrt(2)*upper), moment_residual=residual,
                orthonormal_moment_residual=basis_residual, scaled_basis_condition=condition,
                margins=margins.tolist(), dual=p.tolist(), moment_multipliers=multipliers,
                z_real=z.real.tolist(), z_imag=z.imag.tolist())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository", type=Path)
    parser.add_argument("--nodes", type=int, default=24)
    parser.add_argument("--samples", type=int, default=1024)
    parser.add_argument("--padding", type=int, default=16)
    parser.add_argument("--vertices", type=int, default=256)
    args = parser.parse_args()
    if args.nodes < 8 or args.samples < 256 or args.padding < 8 or args.vertices < 16:
        raise ValueError("Insufficient quadrature or disk resolution")
    parent = load_parent(args.repository)
    prior = parent.load_coherence().prior(args.repository)
    initial = initial_states(prior, args.repository)
    q, matrix, labels, _ = parent.responses(args.repository, args.nodes, args.samples, args.padding)
    w = 4 * np.polynomial.legendre.leggauss(args.nodes)[1] * 3*q*q/8**3
    fine = source_constraints(q, w, initial, False)
    both = source_constraints(q, w, initial, True)
    fine_rows = [i for i, label in enumerate(labels) if label.startswith("fine/wide/number/")]
    groups = {}
    for name, rows in (("fine_wide_full", fine_rows), ("all_24", list(range(len(labels))))):
        selected = matrix[rows]
        group = dict(rows=rows, fine_moments=bounds(selected, fine, args.vertices),
                     imaginary_only=bounds(selected))
        if name == "fine_wide_full":
            group["individual"] = {labels[i].split("/")[-1]: bounds(matrix[[i]], fine, args.vertices)
                                   for i in rows}
        else:
            group["both_mesh_moments"] = bounds(selected, both, args.vertices)
        groups[name] = group
    print(json.dumps(dict(nodes=args.nodes, samples=args.samples, padding=args.padding, vertices=args.vertices,
                          initial=initial, labels=labels, momentum=q.tolist(), groups=groups),
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
