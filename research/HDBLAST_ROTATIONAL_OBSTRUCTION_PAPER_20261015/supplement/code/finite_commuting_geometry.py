#!/usr/bin/env python3
"""Exact finite-matrix checks for the HDBLAST quadrupole reassessment.

Only Sym_0(n)^2 -> so(n), (A,B) |-> [A,B], is tested. No PDE assertion.
Run normally and with ``python -O``. Checks use explicit exceptions, not assert.
Requires SymPy; no external data, floating point, or historical check count.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import sympy as sp


def comm(a, b):
    return a * b - b * a


def simplify_matrix(a):
    return a.applyfunc(sp.cancel)


def is_zero(a):
    return all(sp.cancel(x) == 0 for x in a)


def same(a, b):
    return is_zero(a - b)


def frobenius_squared(a):
    return sp.expand(sum(x * x for x in a))


def pair_norm_squared(a, b):
    return frobenius_squared(a) + frobenius_squared(b)


def derivative(A, B, a, b):
    return comm(a, B) + comm(A, b)


def sym_edge(n, i, j):
    x = sp.zeros(n)
    x[i, j] = x[j, i] = 1
    return x


def skew_edge(n, i, j):
    x = sp.zeros(n)
    x[i, j], x[j, i] = 1, -1
    return x


def edges(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def symmetric_tracefree_basis(n):
    result = []
    for i in range(n - 1):
        h = sp.zeros(n)
        h[i, i], h[n - 1, n - 1] = 1, -1
        result.append(h)
    return result + [sym_edge(n, i, j) for i, j in edges(n)]


def joint_blocks(A, B):
    """Maximal joint eigenspaces; inputs here are simultaneously diagonal."""
    labels = {}
    for i in range(A.rows):
        labels.setdefault((A[i, i], B[i, i]), []).append(i)
    return list(labels.values())


def block_part(a, blocks):
    result = sp.zeros(a.rows)
    for block in blocks:
        for i in block:
            for j in block:
                result[i, j] = a[i, j]
    return result


def block_obstruction(a, b, blocks):
    return comm(block_part(a, blocks), block_part(b, blocks))


def gap_squared(A, B, i, j):
    return sp.expand((A[i, i] - A[j, j]) ** 2
                     + (B[i, i] - B[j, j]) ** 2)


def range_projection(A, B, S):
    result = sp.zeros(A.rows)
    for i, j in edges(A.rows):
        if gap_squared(A, B, i, j) != 0:
            result[i, j], result[j, i] = S[i, j], S[j, i]
    return result


def minimum_norm_inverse(A, B, S):
    """Moore--Penrose inverse in the product Frobenius metric.

    For a singular derivative, this returns the inverse of the range
    projection; a right inverse on the whole skew space exists only at rank N.
    """
    a, b = sp.zeros(A.rows), sp.zeros(A.rows)
    for i, j in edges(A.rows):
        g2 = gap_squared(A, B, i, j)
        if g2 != 0:
            a[i, j] = a[j, i] = (B[j, j] - B[i, i]) * S[i, j] / g2
            b[i, j] = b[j, i] = (A[i, i] - A[j, j]) * S[i, j] / g2
    return a, b


def conjugation_generator(A, B, a, b):
    K = sp.zeros(A.rows)
    for i, j in edges(A.rows):
        da, db = A[j, j] - A[i, i], B[j, j] - B[i, i]
        g2 = da * da + db * db
        if g2 != 0:
            K[i, j] = (da * a[i, j] + db * b[i, j]) / g2
            K[j, i] = -K[i, j]
    return K


def cayley(K, t):
    I = sp.eye(K.rows)
    return simplify_matrix((I + t * K / 2) * (I - t * K / 2).inv())


def tracefree(a):
    return a - sp.trace(a) * sp.eye(a.rows) / a.rows


def admissible(A, B, a, b):
    blocks = joint_blocks(A, B)
    return (is_zero(derivative(A, B, a, b))
            and is_zero(block_obstruction(a, b, blocks)))


class Checks:
    def __init__(self):
        self.items = []

    def check(self, name, condition, details=None):
        if not bool(condition):
            raise RuntimeError(f"Exact check failed: {name}")
        item = {"name": name, "pass": True}
        if details is not None:
            item["details"] = details
        self.items.append(item)

    def rejects(self, name, incorrect_claim):
        self.check(name, not incorrect_claim, {"kind": "wrong-formula control"})


def derivative_matrix(A, B):
    n = A.rows
    basis = symmetric_tracefree_basis(n)
    columns = []
    for h in basis:
        columns.append([derivative(A, B, h, sp.zeros(n))[i, j]
                        for i, j in edges(n)])
    for h in basis:
        columns.append([derivative(A, B, sp.zeros(n), h)[i, j]
                        for i, j in edges(n)])
    return sp.Matrix.hstack(*[sp.Matrix(x) for x in columns])


def make_velocity(A, B):
    """Nontrivial commuting internal blocks plus common orbit directions."""
    n = A.rows
    blocks = joint_blocks(A, B)
    ad, bd, K = sp.zeros(n), sp.zeros(n), sp.zeros(n)
    for block in blocks:
        m = len(block)
        J = sp.zeros(m)
        for i, j in edges(m):
            J[i, j], J[j, i] = i + j + 1, -(i + j + 1)
        U = cayley(J, sp.Rational(1, 2))
        aa = U * sp.diag(*range(1, m + 1)) * U.T
        bb = U * sp.diag(*[k * k for k in range(1, m + 1)]) * U.T
        for i, row in enumerate(block):
            for j, col in enumerate(block):
                ad[row, col], bd[row, col] = aa[i, j], bb[i, j]
    ad, bd = tracefree(ad), tracefree(bd)
    for i, j in edges(n):
        if gap_squared(A, B, i, j) != 0:
            K[i, j], K[j, i] = i + j + 1, -(i + j + 1)
    return ad + comm(K, A), bd + comm(K, B), ad, bd, K


def check_stratum(checks, label, A, B):
    n = A.rows
    blocks = joint_blocks(A, B)
    multiplicities = [len(block) for block in blocks]
    N = n * (n - 1) // 2
    dimV = n * (n + 1) // 2 - 1
    defect = sum(m * (m - 1) // 2 for m in multiplicities)
    rank = int(derivative_matrix(A, B).rank())
    checks.check(f"{label}: commuting tracefree symmetric base",
                 same(A, A.T) and same(B, B.T) and sp.trace(A) == 0
                 and sp.trace(B) == 0 and is_zero(comm(A, B)))
    checks.check(f"{label}: rank and kernel dimension", rank == N - defect,
                 {"joint_multiplicities": multiplicities, "rank": rank,
                  "kernel_dimension": 2 * dimV - rank})

    # Orthonormal off-diagonal input/output bases cancel their common sqrt(2).
    off_columns = []
    for factor in (0, 1):
        for i, j in edges(n):
            h = sym_edge(n, i, j)
            a, b = (h, sp.zeros(n)) if factor == 0 else (sp.zeros(n), h)
            off_columns.append(sp.Matrix([derivative(A, B, a, b)[k, l]
                                          for k, l in edges(n)]))
    Loff = sp.Matrix.hstack(*off_columns)
    gap2s = [gap_squared(A, B, i, j) for i, j in edges(n)]
    checks.check(f"{label}: exact singular-value squares",
                 same(Loff * Loff.T, sp.diag(*gap2s)),
                 {"singular_value_squares": [str(x) for x in gap2s]})

    # Every skew target checks signs, singular range projection, and norm.
    for i, j in edges(n):
        S = skew_edge(n, i, j)
        aa, bb = minimum_norm_inverse(A, B, S)
        checks.check(f"{label}: inverse/projection edge {i + 1}{j + 1}",
                     same(derivative(A, B, aa, bb), range_projection(A, B, S)))
    positive = [g2 for g2 in gap2s if g2 != 0]
    if positive:
        min2 = min(positive)
        i, j = next((i, j) for i, j in edges(n)
                    if gap_squared(A, B, i, j) == min2)
        S = skew_edge(n, i, j)
        aa, bb = minimum_norm_inverse(A, B, S)
        checks.check(f"{label}: sharp inverse norm attained",
                     sp.cancel(pair_norm_squared(aa, bb)
                               / frobenius_squared(S) - 1 / min2) == 0,
                     {"operator_norm_squared": str(1 / min2)})
    else:
        aa, bb = minimum_norm_inverse(A, B, skew_edge(n, 0, 1))
        checks.check(f"{label}: rank-zero pseudoinverse is zero",
                     is_zero(aa) and is_zero(bb))

    a, b, ad, bd, K = make_velocity(A, B)
    checks.check(f"{label}: constructed velocity meets both criteria",
                 admissible(A, B, a, b))
    checks.check(f"{label}: generator reconstructed exactly",
                 same(K, conjugation_generator(A, B, a, b)))
    checks.check(f"{label}: projected full commutator equals internal obstruction",
                 same(block_part(comm(a, b), blocks), block_obstruction(a, b, blocks)))
    c = 2 * comm(K, ad) + comm(K, comm(K, A))
    d = 2 * comm(K, bd) + comm(K, comm(K, B))
    checks.check(f"{label}: explicit analytic-curve second derivative",
                 is_zero(derivative(A, B, c, d) + 2 * comm(a, b)))
    cp, dp = minimum_norm_inverse(A, B, -2 * comm(a, b))
    checks.check(f"{label}: canonical unrestricted second correction",
                 is_zero(derivative(A, B, cp, dp) + 2 * comm(a, b)))

    # Polynomial identities verify the rational curve for all t; no sampling
    # is used to infer an identity. Direct rational evaluations are extra checks.
    t = sp.Symbol("t", real=True)
    M, Nmat = sp.eye(n) - t * K / 2, sp.eye(n) + t * K / 2
    denominator = M.det()
    numerator = Nmat * M.adjugate()
    checks.check(f"{label}: Cayley polynomial orthogonality identity",
                 is_zero((numerator.T * numerator
                          - denominator ** 2 * sp.eye(n)).applyfunc(sp.expand)))
    checks.check(f"{label}: internal polynomial commutator vanishes identically",
                 is_zero(comm(A + t * ad, B + t * bd).applyfunc(sp.expand)))
    Q0 = cayley(K, 0)
    Qprime0 = K  # From (I+tK/2)(I-tK/2)^(-1), checked by polynomial derivative.
    checks.check(f"{label}: Cayley initial value and derivative",
                 same(Q0, sp.eye(n)) and same(
                     numerator.diff(t).subs(t, 0)
                     - numerator.subs(t, 0) * sp.diff(denominator, t).subs(t, 0),
                     Qprime0))
    for value in (sp.Rational(-2), sp.Rational(1, 3), sp.Rational(3, 2)):
        Q = cayley(K, value)
        At, Bt = Q * (A + value * ad) * Q.T, Q * (B + value * bd) * Q.T
        checks.check(f"{label}: exact rational commuting curve t={value}",
                     same(Q.T * Q, sp.eye(n)) and same(At, At.T) and same(Bt, Bt.T)
                     and sp.trace(At) == 0 and sp.trace(Bt) == 0
                     and is_zero(comm(At, Bt)))

    # Equivariance tests ensure the results are not artifacts of diagonal axes.
    J = sum((skew_edge(n, i, j) for i, j in edges(n)), sp.zeros(n))
    U = cayley(J, sp.Rational(1, 3))
    conjugate = lambda x: U * x * U.T
    checks.check(f"{label}: exact orthogonal equivariance",
                 same(derivative(conjugate(A), conjugate(B), conjugate(a), conjugate(b)),
                      conjugate(derivative(A, B, a, b))))


def run_checks():
    checks = Checks()
    strata = [
        ("n4_1111", sp.diag(-3, -1, 1, 3), sp.diag(1, -1, -1, 1)),
        ("n4_211", sp.diag(0, 0, 1, -1), sp.diag(1, 1, -1, -1)),
        ("n4_22", sp.diag(-1, -1, 1, 1), sp.diag(1, 1, -1, -1)),
        ("n4_31", sp.diag(-1, -1, -1, 3), sp.zeros(4)),
        ("n4_4", sp.zeros(4), sp.zeros(4)),
    ]
    for label, A, B in strata:
        check_stratum(checks, label, A, B)

    D = sp.diag(1, -1, 0, 0)
    a, b = sp.diag(0, 0, 1, -1), sym_edge(4, 2, 3)
    checks.check("audit: old seed has rank five", derivative_matrix(D, D).rank() == 5)
    checks.check("audit: both separate singular velocities integrate",
                 admissible(D, D, a, sp.zeros(4)) and admissible(D, D, sp.zeros(4), b))
    checks.check("audit: summed velocity remains in the linearized kernel",
                 is_zero(derivative(D, D, a, b)))
    checks.check("audit: obstruction is exactly 2 in component 34", comm(a, b)[2, 3] == 2)
    checks.rejects("control: linearized kernel alone does not imply a curve",
                   admissible(D, D, a, b))
    t = sp.Symbol("t")
    checks.check("audit: exact straight-path t-squared obstruction",
                 sp.expand(comm(D + t * a, D + t * b)[2, 3]) == 2 * t ** 2)
    checks.check("audit: every linear correction misses component 34",
                 all(derivative(D, D, h, sp.zeros(4))[2, 3] == 0
                     and derivative(D, D, sp.zeros(4), h)[2, 3] == 0
                     for h in symmetric_tracefree_basis(4)))

    R = sp.diag(3, 1, -1, -3) / sp.sqrt(10)
    checks.check("audit: regular replacement has exactly the old norm",
                 sp.trace(D * D) == 2 and sp.trace(R * R) == 2)
    checks.check("audit: regular replacement has rank six", derivative_matrix(R, R).rank() == 6)
    S = skew_edge(4, 0, 1)
    aa, bb = minimum_norm_inverse(R, R, S)
    checks.check("audit: optimal norm-matched inverse constant squared is 5/4",
                 sp.cancel(pair_norm_squared(aa, bb) / frobenius_squared(S)) == sp.Rational(5, 4))
    checks.rejects("control: omitting the Frobenius sqrt-two gives a false bound",
                   pair_norm_squared(aa, bb) <= sp.Rational(5, 8) * frobenius_squared(S))
    h = sym_edge(4, 0, 1)
    checks.check("audit: equal-factor directions vanish at A=B",
                 is_zero(derivative(R, R, h, h)))
    checks.rejects("control: wrong derivative commutator sign",
                   same(derivative(R, R, h, h), comm(h, R) - comm(R, h)))
    checks.rejects("control: wrong inverse sign fails its defining equation",
                   same(derivative(R, R, -aa, bb), S))
    g = sp.sqrt(gap_squared(R, R, 0, 1))
    checks.rejects("control: dividing by joint gap instead of gap squared",
                   same(derivative(R, R, aa * g, bb * g), S))

    A, B = sp.diag(-1, -1, 1, 1), sp.diag(-1, 1, -1, 1)
    checks.check("pencil: both individual spectra repeat yet derivative is onto",
                 len(set(A.diagonal())) == 2 and len(set(B.diagonal())) == 2
                 and derivative_matrix(A, B).rank() == 6)
    checks.check("pencil: explicit generic pencil is simple", len(set((A + 2 * B).diagonal())) == 4)
    checks.rejects("control: repeated A alone does not force singularity",
                   derivative_matrix(A, B).rank() < 6)
    for label, A, B in strata:
        bad_slopes = set()
        for i, j in edges(4):
            db = B[j, j] - B[i, i]
            if db != 0:
                bad_slopes.add(-(A[j, j] - A[i, i]) / db)
        slope = next(k for k in range(len(bad_slopes) + 1) if k not in bad_slopes)
        checks.check(f"{label}: a generic pencil preserves precisely the joint blocks",
                     len(set((A + slope * B).diagonal())) == len(joint_blocks(A, B)))

    A, B = sp.diag(1, -1, 0, 0), sp.zeros(4)
    a, b = sym_edge(4, 0, 1), A
    checks.check("velocity criterion: nonzero full [a,b] can be allowed",
                 admissible(A, B, a, b) and not is_zero(comm(a, b)))
    checks.rejects("control: demanding global [a,b]=0 overrestricts velocities",
                   is_zero(comm(a, b)))
    checks.rejects("control: splitting a maximal repeated joint block is invalid",
                   admissible(sp.zeros(4), sp.zeros(4),
                              sp.diag(1, -1, 0, 0), sym_edge(4, 0, 1)))

    # An arbitrary linearized velocity obeys the off-block cancellation identity;
    # its in-block data need not commute, making second-order solvability fail.
    for label, A, B in strata:
        blocks = joint_blocks(A, B)
        repeated = next((block for block in blocks if len(block) > 1), None)
        if repeated is not None:
            a, b = sp.zeros(4), sp.zeros(4)
            i, j = repeated[:2]
            a[i, i], a[j, j], b[i, j], b[j, i] = 1, -1, 1, 1
            checks.check(f"{label}: nonlinear tangent cone at every repeated stratum",
                         admissible(A, B, a, sp.zeros(4))
                         and admissible(A, B, sp.zeros(4), b)
                         and is_zero(derivative(A, B, a, b))
                         and not admissible(A, B, a, b))
            cp, dp = minimum_norm_inverse(A, B, -2 * comm(a, b))
            checks.rejects(f"{label}: no second correction can reach the obstruction block",
                           is_zero(derivative(A, B, cp, dp) + 2 * comm(a, b)))

    return checks.items


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    items = run_checks()
    receipt = {
        "scope": "Exact finite commuting symmetric-matrix geometry only; no PDE or novelty claim",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "python_optimization": sys.flags.optimize,
        "check_count": len(items),
        "all_pass": all(item["pass"] for item in items),
        "checks": items,
    }
    output = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(json.dumps({key: value for key, value in receipt.items() if key != "checks"}))


if __name__ == "__main__":
    main()
