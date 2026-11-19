"""Bound common separable preparations and evaluate an explicit smooth phase on the archived walls."""
import argparse
import importlib.util
import io
import itertools
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.optimize import brentq


GRIDS = ("fine", "coarse")
WINDOWS = ("wide", "narrow")
BANDS = ("number", "screened_number")
ARMS = ("plus", "minus", "far")


def load_coherence():
    path = Path(__file__).resolve().parent.parent / "coherence/coherence.py"
    spec = importlib.util.spec_from_file_location("phase_response", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def responses(repo, nodes, samples, padding):
    coherence = load_coherence()
    prior = coherence.prior(repo)
    points, quadrature_weights = np.polynomial.legendre.leggauss(nodes)
    q = 4 * (1 + points)
    weights = 4 * quadrature_weights * q**2 / (8**3 / 3)
    matrix, labels, incoherent = [], [], []
    with zipfile.ZipFile(io.BytesIO(prior.checked(repo / prior.WALL, prior.WALL_SHA))) as archive:
        for grid in GRIDS:
            reference = json.loads(prior.checked(repo / prior.REF.format(grid), prior.REF_SHA[grid]))
            scale = reference["matter"]["events"][0]["A_star"]
            results = {}
            for arm in ("quiet", *ARMS):
                geom, start = prior.wall(archive, grid, arm, scale)
                results[arm], _ = coherence.spectrum(geom, start, q, prior.SWITCHES, samples, padding)
            for window, band in itertools.product(WINDOWS, BANDS):
                quiet = results["quiet"][window][band]
                for arm in ARMS:
                    data = results[arm][window][band]
                    difference = weights @ (data[:, 0].real - quiet[:, 0].real)
                    if abs(difference) < 1e-12:
                        raise ValueError("Normalization needs a nonzero incoherent contrast")
                    labels.append(f"{grid}/{window}/{band}/{arm}")
                    incoherent.append(float(difference))
                    matrix.append(-2 * weights * (data[:, 1] - quiet[:, 1]) / difference)
    return q, np.asarray(matrix), labels, np.asarray(incoherent)


def align(vector):
    return np.divide(vector.conj(), abs(vector), out=np.zeros_like(vector), where=abs(vector) > 0)


def bound(matrix):
    """Seek a dual supported on one or two rows, then certify it against every primal row."""
    row_count = len(matrix)
    single = int(np.argmin(abs(matrix).sum(axis=1)))
    p = np.eye(row_count)[single]
    upper = abs(matrix[single]).sum()
    if upper - min((matrix @ align(matrix[single])).real) > 1e-11:
        for first, second in itertools.combinations(range(row_count), 2):
            base = matrix[first]
            delta = matrix[second] - base

            def slope(t):
                return np.real(delta @ align(base + t * delta))

            if slope(0) >= 0 or slope(1) <= 0:
                continue
            mix = brentq(slope, 0., 1., xtol=1e-14)
            candidate = abs(base + mix * delta).sum()
            if candidate < upper:
                upper = candidate
                p = np.zeros(row_count)
                p[first], p[second] = 1 - mix, mix
    combined = p @ matrix
    z = align(combined)
    margins = (matrix @ z).real
    lower, upper = float(margins.min()), float(abs(combined).sum())
    if abs(lower - upper) > 2e-9 or min(p) < 0 or abs(p.sum() - 1) > 1e-14:
        raise ValueError(f"Simplex and state failed to certify this grid: {lower}, {upper}")
    if max(abs(z)) > 1 + 1e-14:
        raise ValueError("Preparation leaves the separable disk")
    return dict(lower=lower, upper=upper, margin_factors=margins.tolist(), dual=p.tolist(),
                z_real=z.real.tolist(), z_imag=z.imag.tolist())


def common_angle_row_ceiling(vector):
    """Farthest vertex of the planar sum of line segments [0,vector[q]]."""
    vectors = vector[abs(vector) > 0]
    angles = np.angle(vectors)
    boundaries = np.unique(np.mod(np.r_[np.pi/2 - angles, -np.pi/2 - angles], 2*np.pi))
    ends = np.r_[boundaries[1:], boundaries[0] + 2*np.pi]
    probes = (boundaries + ends) / 2
    values = []
    for theta in probes:
        polar = np.exp(1j * theta)
        active = (vectors * polar).real > 0
        values.append(abs(vectors[active].sum()))
    return float(max(values))


def subsets(labels):
    keys = sorted(set(label.rsplit("/", 1)[0] for label in labels))
    groups = {key: [i for i, label in enumerate(labels) if label.startswith(key + "/")] for key in keys}
    for grid in GRIDS:
        groups[grid + "/all"] = [i for i, label in enumerate(labels) if label.startswith(grid + "/")]
    groups["all"] = list(range(len(labels)))
    return groups


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository", type=Path)
    parser.add_argument("--nodes", type=int, default=24)
    parser.add_argument("--samples", type=int, default=1024)
    parser.add_argument("--padding", type=int, default=16)
    args = parser.parse_args()
    if args.nodes < 4 or args.samples < 256 or args.padding < 8:
        raise ValueError("Insufficient response sampling")
    q, matrix, labels, incoherent = responses(args.repository, args.nodes, args.samples, args.padding)
    groups = {}
    for name, rows in subsets(labels).items():
        selected = matrix[rows]
        groups[name] = dict(rows=rows, free=bound(selected),
                            common=bound(selected.sum(axis=1, keepdims=True)))

    z = np.exp(1j * (q*q/20 - 2.43))
    margins = (matrix @ z).real
    smooth = dict(phase="q^2/20 - 2.43", minimum=float(margins.min()), maximum=float(margins.max()),
                  margin_factors=margins.tolist(), contrasts=(incoherent * (1 - margins)).tolist())
    ceilings = [common_angle_row_ceiling(row) for row in matrix]
    limiting = int(np.argmin(ceilings))
    common_angle = dict(upper=ceilings[limiting], row=limiting, row_upper_bounds=ceilings)
    out = dict(momentum=q.tolist(), incoherent=incoherent.tolist(), labels=labels, nodes=args.nodes,
               samples=args.samples, padding=args.padding, groups=groups, smooth=smooth,
               common_angle_variable_modulus=common_angle)
    print(json.dumps(out, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
