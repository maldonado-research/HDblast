#!/usr/bin/env python3
"""All-point sphere-curl audit in the frozen selected tensor basis.

Standard library only.  Primary quotient arithmetic is frozen input; the
independent route uses unreduced ambient polynomials, three-slot tangent
projection, and homogenization, with a separate sphere moment formula.
"""
import argparse
import importlib.util
import json
import os
import sys
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations, product
from math import factorial
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
INPUT_SHA256 = "2c0bafec777581c79c3fd6bd9879d0d3cb8114bb899722bd33882086751ded60"
I = range(4)
Z = (0, 0, 0, 0)
ONE = {Z: Q(1)}
CONTROLS = {
    "orientation_flip": "CURL_MAP",
    "symmetrization_factor": "CURL_MAP",
    "drop_minus": "CURL_MAP",
    "equal_cancellation": "EQUAL_COMBINATION",
    "projection_factor": "TT_PROJECTION_NORMALIZATION",
    "independent_epsilon_flip": "INDEPENDENT_CURL_MAP",
    "quadratic_factor": "GRADIENT_HESSIAN_IDENTITY",
}


def need(ok, label):
    if not ok:
        raise ValueError(label)


def fresh_output_path(value):
    """Reject unsafe destinations before any project input or calculation.

    Final exclusive creation also prevents following a newly introduced final
    symlink.  Existing parent directories are required; symlink components
    are rejected so the reviewed lexical destination matches its target.
    """
    path = Path(os.path.abspath(value))
    for part in (path,) + tuple(path.parents):
        need(not part.is_symlink(), "OUTPUT_SYMLINK")
    need(path.parent.is_dir(), "OUTPUT_PARENT_REQUIRED")
    resolved = path.resolve()
    protected = [HERE.resolve()]
    protected.extend(p.resolve() for p in (HERE,) + tuple(HERE.parents)
                     if p.name == "HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026")
    need(all(resolved != p and p not in resolved.parents for p in protected),
         "OUTPUT_PROTECTED_SOURCE")
    # lexists, unlike exists, also sees dangling symlinks; no file is opened.
    need(not os.path.lexists(path), "OUTPUT_EXISTS")
    return path


def load_inputs(input_path=None):
    path = Path(input_path) if input_path else HERE / "inputs" / "angular_exact.py"
    need(sha256(path.read_bytes()).hexdigest() == INPUT_SHA256, "FROZEN_INPUT_PIN")
    spec = importlib.util.spec_from_file_location("screen_frozen_angular", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def curl_tensor(t, a=None, orientation=1):
    """Return epsilon_(i^{cd} D_c t_{j)d}; normalized symmetrization.

    ``t`` has frozen ambient polynomial tangent slots.  Orientation +1 means
    epsilon_ijk=n^A EPS_Aijk with EPS_0123=+1.  This API does not project its
    result onto a presumed finite harmonic feature space.
    """
    a = load_inputs() if a is None else a
    dt = a.covd(t)
    eps = {(i, c, d): a.add(*(a.scale(orientation * a.EPS.get((b, i, c, d), 0), a.X[b])
                            for b in I)) for i, c, d in product(I, repeat=3)}
    return {(i, j): a.scale(Q(1, 2), a.add(
        *(a.mul(eps[i, c, d], dt[c, j, d]) for c, d in product(I, repeat=2)),
        *(a.mul(eps[j, c, d], dt[c, i, d]) for c, d in product(I, repeat=2))))
        for i, j in product(I, repeat=2)}


def project_tt_pair(t, a=None):
    """Return the exact two selected amplitudes 6 <t,T+>, 6 <t,T->.

    This is a pair projection, not a claim that arbitrary t has no additional
    TT harmonics.  No reconstruction hypothesis is required for these two
    inner products.
    """
    a = load_inputs() if a is None else a
    return tuple(6 * a.inner(t, b) for b in a.TT)


def gradient_hessian_curl(a=None):
    """epsilon_(i^{cd} (D_cY) (D_j)D_dY) in frozen ambient slots."""
    a = load_inputs() if a is None else a
    grad, hess = a.grad(a.Y), a.hess(a.Y)
    eps = {(i, c, d): a.add(*(a.scale(a.EPS.get((b, i, c, d), 0), a.X[b])
                            for b in I)) for i, c, d in product(I, repeat=3)}
    return {(i, j): a.scale(Q(1, 2), a.add(
        *(a.mul(eps[i, c, d], a.mul(grad[c,], hess[j, d])) for c, d in product(I, repeat=2)),
        *(a.mul(eps[j, c, d], a.mul(grad[c,], hess[i, d])) for c, d in product(I, repeat=2))))
        for i, j in product(I, repeat=2)}


# Independent raw polynomial arithmetic: no sphere reduction at any stage.
def radd(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return {e: c for e, c in out.items() if c}


def rscale(c, p):
    return {e: c * v for e, v in p.items() if c * v}


def rmul(*ps):
    out = ONE
    for p in ps:
        nxt = {}
        for e, c in out.items():
            for f, d in p.items():
                g = tuple(x + y for x, y in zip(e, f))
                nxt[g] = nxt.get(g, Q(0)) + c * d
        out = {e: c for e, c in nxt.items() if c}
    return out


def rdiff(p, k):
    out = {}
    for e, c in p.items():
        if e[k]:
            f = list(e)
            f[k] -= 1
            out[tuple(f)] = c * e[k]
    return out


RX = [{tuple(int(j == i) for j in I): Q(1)} for i in I]
RADIUS2 = radd(*(rmul(x, x) for x in RX))
RP = {(i, j): radd(ONE if i == j else {}, rscale(-1, rmul(RX[i], RX[j])))
      for i, j in product(I, repeat=2)}


def parity(p):
    return (-1) ** sum(p[j] > p[k] for j in I for k in range(j + 1, 4))


REPS = {p: parity(p) for p in permutations(I)}


def homogenized(p):
    """Exact zero is sufficient for p|S3=0; no points are evaluated.

    Every polynomial checked here has one degree parity. Homogenize its
    terms to the maximum degree by multiplying by powers of radius squared.
    Since radius squared=1 on S3, a zero homogenization proves the identity.
    """
    if not p:
        return {}
    degree = max(map(sum, p))
    need(all((degree - sum(e)) % 2 == 0 for e in p), "HOMOGENIZATION_PARITY")
    powers = [ONE]
    for _ in range(degree // 2):
        powers.append(rmul(powers[-1], RADIUS2))
    return radd(*(rmul({e: c}, powers[(degree - sum(e)) // 2]) for e, c in p.items()))


def raw_expect(t, u, label):
    need(all(not homogenized(radd(t[k], rscale(-1, u[k]))) for k in t), label)


def raw_frozen_tt():
    # An independent bivector-diagonal definition, rather than a.C or a.make_tt.
    section = {(0, 1): Q(-2, 3), (2, 3): Q(-2, 3), (0, 2): Q(1, 3),
               (0, 3): Q(1, 3), (1, 2): Q(1, 3), (1, 3): Q(1, 3)}

    def C(i, j, k, l):
        if i == j or k == l or sorted((i, j)) != sorted((k, l)):
            return Q(0)
        return section[tuple(sorted((i, j)))] * (1 if (i < j) == (k < l) else -1)

    out = []
    for s in (1, -1):
        t = {}
        for i, j in product(I, repeat=2):
            terms = []
            for b, d in product(I, repeat=2):
                star = sum((REPS.get((i, b, c, e), 0) * C(c, e, j, d)
                            for c, e in product(I, repeat=2)), Q(0)) / 2
                terms.append(rscale((C(i, b, j, d) + s * star) / 2, rmul(RX[b], RX[d])))
            t[i, j] = radd(*terms)
        out.append(t)
    return out


def raw_curl(t, orientation=1):
    # Fully projected ambient ordinary derivative, independently of covd's
    # connection formula or normal-form polynomial representatives.
    ordinary = {(b, i, j): rdiff(t[i, j], b) for b, i, j in product(I, repeat=3)}
    dt = {}
    for c, j, d in product(I, repeat=3):
        dt[c, j, d] = radd(*(rmul(RP[c, b], RP[j, e], RP[d, f], ordinary[b, e, f])
                            for b, e, f in product(I, repeat=3)))
    eps = {(i, c, d): radd(*(rscale(orientation * REPS.get((b, i, c, d), 0), RX[b])
                            for b in I)) for i, c, d in product(I, repeat=3)}
    return {(i, j): rscale(Q(1, 2), radd(
        *(rmul(eps[i, c, d], dt[c, j, d]) for c, d in product(I, repeat=2)),
        *(rmul(eps[j, c, d], dt[c, i, d]) for c, d in product(I, repeat=2))))
        for i, j in product(I, repeat=2)}


def raw_dot(t, u):
    return radd(*(rmul(t[k], u[k]) for k in t))


def raw_gradient_hessian():
    y = radd(rmul(RX[0], RX[0]), rscale(-1, rmul(RX[1], RX[1])))
    grad = {i: radd(*(rmul(RP[i, b], rdiff(y, b)) for b in I)) for i in I}
    hess = {(i, j): radd(*(rmul(RP[i, b], RP[j, d], rdiff(grad[d], b))
                           for b, d in product(I, repeat=2))) for i, j in product(I, repeat=2)}
    eps = {(i, c, d): radd(*(rscale(REPS.get((b, i, c, d), 0), RX[b])
                            for b in I)) for i, c, d in product(I, repeat=3)}
    return {(i, j): rscale(Q(1, 2), radd(
        *(rmul(eps[i, c, d], grad[c], hess[j, d]) for c, d in product(I, repeat=2)),
        *(rmul(eps[j, c, d], grad[c], hess[i, d]) for c, d in product(I, repeat=2))))
        for i, j in product(I, repeat=2)}


def raw_mean(p):
    # E[n_0^(2k0)...n_3^(2k3)] = prod((2ki)!/(4^ki ki!))/(m+1)!.
    total = Q(0)
    for e, c in p.items():
        if any(k % 2 for k in e):
            continue
        ks = [k // 2 for k in e]
        value = Q(1, factorial(sum(ks) + 1))
        for k in ks:
            value *= Q(factorial(2 * k), 4 ** k * factorial(k))
        total += c * value
    return total


def maximum_degree(t):
    return max((sum(e) for p in t.values() for e in p), default=0)


def run(inject=None, input_path=None):
    a = load_inputs(input_path)
    zero = a.zero_tensor(2)
    expected = [zero] * 5 + [a.tscale(3, a.TT[0]), a.tscale(-3, a.TT[1])]
    derived = [curl_tensor(t, a, -1 if inject == "orientation_flip" else 1) for t in a.TENSORS]
    if inject == "symmetrization_factor":
        derived = [a.tscale(2, t) for t in derived]
    if inject == "drop_minus":
        derived[-1] = zero
    for t, u in zip(derived, expected):
        a.expect(t, u, "CURL_MAP")

    # The frozen seven tensor columns have symmetric, tangential outputs.
    for t in derived:
        a.expect(t, {(i, j): t[j, i] for i, j in t}, "CURL_SYMMETRY")
        a.expect(a.trace(t), a.zero_tensor(0), "CURL_TRACE")
        a.expect(a.div(t), a.zero_tensor(1), "CURL_DIVERGENCE")
        for j in I:
            need(not a.add(*(a.mul(a.X[i], t[i, j]) for i in I)), "CURL_TANGENCY")

    raw_tt = raw_frozen_tt()
    for t, frozen in zip(raw_tt, a.TT):
        raw_expect(t, frozen, "INDEPENDENT_FROZEN_BASIS")
    raw_basis = list(a.TENSORS[:5]) + raw_tt
    raw_derived = [raw_curl(t, -1 if inject == "independent_epsilon_flip" else 1) for t in raw_basis]
    for t, u in zip(raw_derived, expected):
        raw_expect(t, u, "INDEPENDENT_CURL_MAP")

    # All norm and projection constants are separately integrated in raw form.
    gram = [[raw_mean(raw_dot(t, u)) for u in raw_tt] for t in raw_tt]
    need(gram == [[Q(1, 6), 0], [0, Q(1, 6)]], "INDEPENDENT_GRAM")
    for i, t in enumerate(a.TT):
        coeff = list(project_tt_pair(t, a))
        if inject == "projection_factor":
            coeff = [c / 2 for c in coeff]
        need(coeff == [Q(i == j) for j in range(2)], "TT_PROJECTION_NORMALIZATION")
    for t in a.TENSORS[:5]:
        need(project_tt_pair(t, a) == (0, 0), "SCALAR_TT_ORTHOGONALITY")

    plus, minus = a.TT
    summed = a.tadd(plus, minus)
    difference = a.tadd(plus, a.tscale(-1, minus))
    equal_curl = curl_tensor(summed, a)
    if inject == "equal_cancellation":
        equal_curl = zero
    a.expect(equal_curl, a.tscale(3, difference), "EQUAL_COMBINATION")
    S = a.S
    norm = a.scale(18, a.mul(S, a.add(a.ONE, a.scale(-1, S))))
    need(a.dot(equal_curl, equal_curl) == norm, "POINTWISE_EQUAL_CURL_NORM")
    need(a.inner(equal_curl, equal_curl) == 3, "INTEGRATED_EQUAL_CURL_NORM")

    rS = radd(rmul(RX[0], RX[0]), rmul(RX[1], RX[1]))
    for t in raw_tt:
        need(not homogenized(radd(raw_dot(t, t), rscale(-Q(1, 6), ONE))), "INDEPENDENT_POINTWISE_TT_NORM")
    cross = radd(rscale(Q(1, 6), ONE), rscale(-1, rS), rmul(rS, rS))
    need(not homogenized(radd(raw_dot(*raw_tt), rscale(-1, cross))), "INDEPENDENT_POINTWISE_CROSS")
    raw_diff = {k: radd(raw_tt[0][k], rscale(-1, raw_tt[1][k])) for k in raw_tt[0]}
    raw_norm = rscale(18, rmul(rS, radd(ONE, rscale(-1, rS))))
    need(not homogenized(radd(rscale(9, raw_dot(raw_diff, raw_diff)), rscale(-1, raw_norm))),
         "INDEPENDENT_EQUAL_CURL_NORM")
    need(raw_mean(raw_norm) == 3, "INDEPENDENT_INTEGRATED_NORM")

    gh = gradient_hessian_curl(a)
    if inject == "quadratic_factor":
        gh = a.tscale(Q(1, 2), gh)
    a.expect(gh, a.tscale(2, difference), "GRADIENT_HESSIAN_IDENTITY")
    raw_gh = raw_gradient_hessian()
    raw_expect(raw_gh, a.tscale(2, difference), "INDEPENDENT_GRADIENT_HESSIAN_IDENTITY")
    need(a.inner(gh, gh) == Q(4, 3), "GRADIENT_HESSIAN_NORM")

    matrix = [[Q(0)] * 7 for _ in range(7)]
    matrix[5][5] = Q(3)
    matrix[6][6] = Q(-3)
    return {
        "status": "PASS_EXACT_FROZEN_SPHERE_CURL",
        "orientation": "EPS_0123=+1; epsilon_ijk=n^A EPS_Aijk; at n=e0, epsilon_123=+1",
        "tensor_basis": a.TENSOR_NAMES,
        "curl_matrix": [[str(c) for c in row] for row in matrix],
        "TT_curl_matrix": [["3", "0"], ["0", "-3"]],
        "TT_Gram_probability_mean": [[str(c) for c in row] for row in gram],
        "TT_projection_amplitudes": ["6*mean(H:Tplus)", "6*mean(H:Tminus)"],
        "pointwise": {
            "S": "n0^2+n1^2", "norm_Tplus_squared": "1/6", "norm_Tminus_squared": "1/6",
            "Tplus_colon_Tminus": "1/6-S*(1-S)",
            "curl_equal_combination": "3*(Tplus-Tminus)",
            "curl_equal_norm_squared": "18*S*(1-S)",
            "curl_equal_zero_set": "S=0 or S=1 (two coordinate great circles)",
        },
        "general_amplitudes": {
            "curl": "3*A*Tplus-3*B*Tminus",
            "curl_mean_norm_squared": "(3/2)*(A^2+B^2)",
            "curl_pointwise_norm_squared": "(3/2)*(A-B)^2+18*A*B*S*(1-S)",
        },
        "quadratic_angular_identity": {
            "definition": "A_ab=epsilon_(a^{cd} (D_cY) (D_b)D_dY)",
            "Y": "n0^2-n1^2",
            "complete_tensor": "A=2*Tplus-2*Tminus",
            "TT_pair_amplitudes": ["2", "-2"],
            "pointwise_norm_squared": "8*S*(1-S)",
            "mean_norm_squared": "4/3",
            "proofs": "quotient zero coefficients and independent raw homogeneous zero coefficients",
        },
        "linear_screen_magnetic_TT": {
            "amplitudes": ["3*uplus_x/(2*ell*rho^2*cosh(t))", "-3*uminus_x/(2*ell*rho^2*cosh(t))"],
            "equal_amplitude_mean_norm_squared": "3*ux^2/(4*ell^2*rho^4*cosh(t)^2)",
            "scope": "linear second-metric term only; quadratic curvature/observer terms separate",
        },
        "all_point_proofs": {
            "primary": "all 7 tensor columns; exact zero polynomials in Q[n]/(sum(n^2)-1)",
            "independent": "raw unreduced polynomials; three-slot tangent derivative; zero homogeneous coefficients",
            "raw_curl_max_degrees": [maximum_degree(t) for t in raw_derived],
            "angular_points_sampled": 0, "ordinary_floating_point_decisions": 0,
            "symmetry_tangency_trace_divergence": True,
        },
        "scope": "Exact selected seven-column curl map and selected TT pair projection, not arbitrary full TT closure or complete H2.",
        "input_sha256": {"angular_exact.py": INPUT_SHA256},
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--input", type=Path, help="pinned frozen angular_exact.py override")
    parser.add_argument("--inject", choices=sorted(CONTROLS))
    args = parser.parse_args()
    output = fresh_output_path(args.output) if args.output else None
    result = run(args.inject, args.input)
    need(args.inject is None, "CONTROL_NOT_REJECTED")
    data = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if output:
        with output.open("x", encoding="utf-8") as stream:
            stream.write(data)
        print(result["status"])
    else:
        print(data, end="")


if __name__ == "__main__":
    main()
