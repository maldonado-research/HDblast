#!/usr/bin/env python3
"""Exact, unpublished S3 evaluation bound for the selected 23 angular features.

Standard library only; all certificate checks are explicit (also under -O).
No evolution/source closure or reduced numerical response is certified here.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations, product
from math import comb
from pathlib import Path
import json

I = range(4)
Z4 = (0, 0, 0, 0)
Z2 = (0, 0)
CONTROLS = {
    "wrong_F_norm": "FULL_ANGULAR_GRAM",
    "omit_mixed_weight": "FROBENIUS_FIELD_WEIGHTS",
    "omit_trace_weight": "FROBENIUS_FIELD_WEIGHTS",
    "ambient_hessian": "SPHERE_TANGENT_TF",
    "drop_TT": "FULL_ANGULAR_GRAM",
    "wrong_moment": "FULL_ANGULAR_GRAM",
    "constant_kernel23": "NONCONSTANT_FULL_KERNEL",
    "unsafe_constant4": "SCALAR_SHARP_ENDPOINT",
    "corrupt_Bernstein": "BERNSTEIN_RECONSTRUCTION",
}


def need(condition, code):
    if not condition:
        raise ValueError(code)


def clean(p):
    return {e: Q(v) for e, v in p.items() if v}


def add(*polys):
    out = {}
    for p in polys:
        for e, v in p.items():
            out[e] = out.get(e, Q(0)) + v
    return clean(out)


def scale(c, p):
    return clean({e: c * v for e, v in p.items()})


def mul(p, q):
    out = {}
    for e, v in p.items():
        for f, w in q.items():
            h = tuple(a + b for a, b in zip(e, f))
            out[h] = out.get(h, Q(0)) + v * w
    return clean(out)


def power(p, n):
    out = {tuple(0 for _ in next(iter(p))): Q(1)}
    for _ in range(n):
        out = mul(out, p)
    return out


def derivative(p, i):
    out = {}
    for e, v in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = v * e[i]
    return clean(out)


def moment(e, mutation=None):
    if any(k % 2 for k in e):
        return Q(0)
    out = Q(1)
    degree = sum(e) // 2
    for k in e:
        for j in range(k // 2):
            out *= 2 * j + 1
    for j in range(degree):
        out /= 4 + 2 * j
    if mutation == "wrong_moment" and degree == 2:
        out *= 2
    return out


def integral(p, mutation=None):
    return sum((v * moment(e, mutation) for e, v in p.items()), Q(0))


def inner(p, q, mutation=None):
    return integral(mul(p, q), mutation)


def vector_inner(u, v, mutation=None):
    return sum((inner(a, b, mutation) for a, b in zip(u, v)), Q(0))


def norm_polynomial(v):
    return add(*(mul(a, a) for a in v))


def sphere_reduce(p):
    """Exact elimination of n4^2; require invariance in the 34-plane.

    Every input is checked even in each coordinate. Put X=n1^2,Y=n2^2,
    W=n3^2 and n4^2=1-X-Y-W, expand, and check W cancels identically.
    The output is thus an exact polynomial on X,Y>=0, X+Y<=1.
    """
    need(all(all(k % 2 == 0 for k in e) for e in p), "EVEN_NORM_POLYNOMIAL")
    one = {(0, 0, 0): Q(1)}
    replacement = {(0, 0, 0): Q(1), (1, 0, 0): Q(-1),
                   (0, 1, 0): Q(-1), (0, 0, 1): Q(-1)}
    out = {}
    for e, value in p.items():
        a, b, c, d = (k // 2 for k in e)
        term = mul({(a, b, c): value}, power(replacement, d) if d else one)
        out = add(out, term)
    need(all(e[2] == 0 for e in out), "CONTINUOUS_34_INVARIANCE")
    return clean({(a, b): value for (a, b, c), value in out.items()})


def bernstein_coefficients(p, degree, mutation=None):
    """Triangle Bernstein coefficients in X,Y,R, X+Y+R=1.

    X^a Y^b has coefficient (i)_a (j)_b / (degree)_(a+b)
    at the degree-degree Bernstein basis indexed by i+j+k=degree.
    """
    need(degree >= max((sum(e) for e in p), default=0), "BERNSTEIN_DEGREE")
    out = {}
    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            value = Q(0)
            for (a, b), c in p.items():
                if i >= a and j >= b:
                    value += c * Q(comb(i, a) * comb(j, b),
                                   comb(degree, a + b) * comb(a + b, a))
            out[i, j, degree - i - j] = value
    if mutation == "corrupt_Bernstein":
        out[0, 0, degree] += 1
    return out


def reconstruct_bernstein(coeffs, degree):
    """Reconstruct in X,Y using R=1-X-Y; no sampling enters the proof."""
    r = {Z2: Q(1), (1, 0): Q(-1), (0, 1): Q(-1)}
    out = {}
    for (i, j, k), coefficient in coeffs.items():
        term = mul({(i, j): coefficient * comb(degree, i) * comb(degree - i, j)},
                   power(r, k) if k else {Z2: Q(1)})
        out = add(out, term)
    return out


def polynomial_json(p):
    return [{"powers_X_Y": list(e), "coefficient": str(v)} for e, v in sorted(p.items())]


def sign(perm):
    return -1 if sum(perm[i] > perm[j] for i in I for j in range(i + 1, 4)) % 2 else 1


def build_features(mutation):
    one = {Z4: Q(1)}
    x = []
    for i in I:
        e = [0] * 4
        e[i] = 1
        x.append({tuple(e): Q(1)})
    xx = [mul(v, v) for v in x]
    radius = add(*xx)
    s = add(xx[0], xx[1])
    y = add(xx[0], scale(-1, xx[1]))
    z = add(s, scale(-Q(1, 2), radius))
    f = add(mul(y, y), scale(-Q(1, 2), mul(radius, s)), scale(Q(1, 12), mul(radius, radius)))
    p = {(i, j): add(one if i == j else {}, scale(-1, mul(x[i], x[j])))
         for i, j in product(I, repeat=2)}
    gradients = []
    hessians = []
    for feature, degree in [(z, 2), (f, 4)]:
        df = [derivative(feature, i) for i in I]
        gradients.append([add(df[i], scale(-degree, mul(feature, x[i]))) for i in I])
        dd = {(i, j): derivative(df[i], j) for i, j in product(I, repeat=2)}
        if mutation == "ambient_hessian":
            hessians.append([dd[i, j] for i, j in product(I, repeat=2)])
        else:
            hessians.append([
                add(*(mul(mul(p[i, a], dd[a, b]), p[b, j]) for a, b in product(I, repeat=2)),
                    scale(Q(degree * (degree - 1), 3), mul(feature, p[i, j])))
                for i, j in product(I, repeat=2)])

    delta = lambda i, j: Q(i == j)
    eigen = [Q(1), Q(-1), Q(0), Q(0)]
    a = lambda i, j: eigen[i] * delta(i, j)
    raw = {(i, j, k, l): 2 * (a(i, k) * a(j, l) - a(i, l) * a(j, k))
           for i, j, k, l in product(I, repeat=4)}
    ricci = {(j, l): sum(raw[i, j, i, l] for i in I) for j, l in product(I, repeat=2)}
    scalar = sum(ricci[i, i] for i in I)
    weyl = {}
    for i, j, k, l in product(I, repeat=4):
        traces = (delta(i, k) * ricci[j, l] - delta(i, l) * ricci[j, k]
                  - delta(j, k) * ricci[i, l] + delta(j, l) * ricci[i, k]) / 2
        weyl[i, j, k, l] = raw[i, j, k, l] - traces + scalar * (delta(i, k) * delta(j, l) - delta(i, l) * delta(j, k)) / 6
    epsilon = {perm: sign(perm) for perm in permutations(I)}
    star = {(i, j, k, l): sum(Q(epsilon.get((i, j, a, b), 0)) * weyl[a, b, k, l] / 2
                                 for a, b in product(I, repeat=2))
            for i, j, k, l in product(I, repeat=4)}
    tt = []
    for chirality in [1, -1]:
        c = {key: (weyl[key] + chirality * star[key]) / 2 for key in weyl}
        tensor = [add(*(scale(c[i, a, j, b], mul(x[a], x[b]))
                        for a, b in product(I, repeat=2))) for i, j in product(I, repeat=2)]
        if mutation == "drop_TT" and chirality == -1:
            tensor = [{} for _ in tensor]
        tt.append(tensor)

    # These identities are polynomial identities on the sphere, checked by
    # direct exact moments of squared residuals. Nonnegative residual norm has
    # zero integral iff the polynomial vanishes on the sphere.
    for tensor in hessians + tt:
        tangent_residual = [add(*(mul(x[i], tensor[4 * i + j]) for i in I)) for j in I]
        trace_residual = add(*(tensor[4 * i + i] for i in I))
        symmetry = [add(tensor[4 * i + j], scale(-1, tensor[4 * j + i]))
                    for i, j in product(I, repeat=2)]
        need(integral(norm_polynomial(tangent_residual + [trace_residual] + symmetry)) == 0,
             "SPHERE_TANGENT_TF")
    # Physical symmetric-tensor slots are y,0 followed by four ambient angular
    # components (restricted tangentially); the final scalar component is w.
    # Constructing the embeddings makes factors 2 and 3 independently visible.
    physical_scalars = []
    for slot in ["yy", "y0", "00", "trace", "w"]:
        embedding = [{} for _ in range(37)]
        if slot == "yy":
            embedding[0] = one
        elif slot == "y0":
            embedding[1] = one
            embedding[6] = one
        elif slot == "00":
            embedding[7] = one
        elif slot == "trace":
            for i, j in product(I, repeat=2):
                embedding[6 * (i + 2) + j + 2] = p[i, j]
        else:
            embedding[36] = one
        physical_scalars.append(embedding)
    geometric_weights = [integral(norm_polynomial(v)) for v in physical_scalars]
    need([sphere_reduce(norm_polynomial(v)) for v in physical_scalars]
         == [{Z2: q} for q in [Q(1), Q(2), Q(1), Q(3), Q(1)]],
         "PHYSICAL_POINTWISE_FIELD_WEIGHTS")
    for i, u in enumerate(physical_scalars):
        for j, v in enumerate(physical_scalars):
            need(vector_inner(u, v) == (geometric_weights[i] if i == j else 0),
                 "PHYSICAL_SCALAR_SLOT_ORTHOGONALITY")
    return one, z, f, gradients, hessians, tt, geometric_weights


def run(mutation=None):
    one, z, f, gradients, hessians, tt, geometric_weights = build_features(mutation)
    scalar_norms = [Q(1), Q(1, 12), Q(19, 720)]
    if mutation == "wrong_F_norm":
        scalar_norms[2] = Q(1, 40)
    gradient_norms = [Q(2, 3), Q(19, 30)]
    tensor_norms = [Q(20, 9), Q(133, 15), Q(1, 6), Q(1, 6)]
    scalar_weights = [Q(1), Q(2), Q(1), Q(3), Q(1)]
    gradient_weights = [Q(2), Q(2)]
    if mutation == "omit_mixed_weight":
        scalar_weights[1] = Q(1)
    if mutation == "omit_trace_weight":
        scalar_weights[3] = Q(1)
    need(scalar_weights == geometric_weights and gradient_weights == [2, 2], "FROBENIUS_FIELD_WEIGHTS")

    scalars = [one, z, f]
    tensors = hessians + tt
    gram_blocks = [(scalars, scalar_norms, lambda a, b: inner(a, b, mutation)),
                   (gradients, gradient_norms, lambda a, b: vector_inner(a, b, mutation)),
                   (tensors, tensor_norms, lambda a, b: vector_inner(a, b, mutation))]
    gram_checks = 0
    for basis, norms, pairing in gram_blocks:
        for i, u in enumerate(basis):
            for j, v in enumerate(basis):
                need(pairing(u, v) == (norms[i] if i == j else 0), "FULL_ANGULAR_GRAM")
                gram_checks += 1
    need(5 * len(scalars) + 2 * len(gradients) + len(tensors) == 23, "FEATURE_DIMENSION")
    # The physical embedding norm and its L2 Gram both multiply by the
    # corresponding Frobenius weight. Thus the normalized kernel is independent
    # of that weight, but only after verifying the correct weights above.
    k_scalar = sphere_reduce(add(*(scale(1 / q, mul(v, v)) for v, q in zip(scalars, scalar_norms))))
    k_gradient = sphere_reduce(add(*(scale(1 / q, norm_polynomial(v)) for v, q in zip(gradients, gradient_norms))))
    k_tensor = sphere_reduce(add(*(scale(1 / q, norm_polynomial(v)) for v, q in zip(tensors, tensor_norms))))
    tt_pointwise = [sphere_reduce(norm_polynomial(v)) for v in tt]
    need(tt_pointwise == [{Z2: Q(1, 6)}, {Z2: Q(1, 6)}], "TT_CONSTANT_POINTWISE_NORM")
    traces = {"scalar_block": k_scalar, "gradient_block": k_gradient, "TF_tensor_block": k_tensor}
    expected_integrals = [Q(3), Q(2), Q(4)]
    for polynomial, expected in zip(traces.values(), expected_integrals):
        lifted = {(2 * i, 2 * j, 0, 0): v for (i, j), v in polynomial.items()}
        need(integral(lifted) == expected, "BLOCK_TRACE_INTEGRALS")

    squared_bound = Q(321, 19)
    if mutation == "unsafe_constant4":
        squared_bound = Q(16)
    at_e1 = sum((value for (i, j), value in k_scalar.items() if j == 0), Q(0))
    need(at_e1 == Q(321, 19) and at_e1 <= squared_bound, "SCALAR_SHARP_ENDPOINT")
    certificates = {}
    # The first successful degree is deterministic and the exact reconstruction
    # is checked independently from the coefficient conversion formula.
    for name, polynomial in traces.items():
        for degree in range(4, 25):
            coefficients = bernstein_coefficients(polynomial, degree, mutation)
            need(reconstruct_bernstein(coefficients, degree) == polynomial, "BERNSTEIN_RECONSTRUCTION")
            upper = max(coefficients.values())
            if upper <= squared_bound:
                certificates[name] = {
                    "degree": degree, "coefficient_count": len(coefficients),
                    "minimum_coefficient": str(min(coefficients.values())),
                    "maximum_coefficient": str(upper),
                    "coefficients": [{"i_j_k": list(index), "value": str(value)}
                                     for index, value in sorted(coefficients.items())],
                    "polynomial": polynomial_json(polynomial),
                }
                break
        need(name in certificates, "GLOBAL_BERNSTEIN_BOUND")
    full_trace = add(scale(5, k_scalar), scale(2, k_gradient), k_tensor)
    full_trace_at_e1 = sum((value for (i, j), value in full_trace.items() if j == 0), Q(0))
    full_trace_at_e3 = full_trace.get(Z2, Q(0))
    if mutation == "constant_kernel23":
        need(full_trace_at_e1 == 23 and full_trace_at_e3 == 23, "NONCONSTANT_FULL_KERNEL")
    need(full_trace_at_e1 != full_trace_at_e3, "NONCONSTANT_FULL_KERNEL")
    need(integral({(2 * i, 2 * j, 0, 0): v for (i, j), v in full_trace.items()}) == 23,
         "FULL_TRACE_INTEGRAL")
    return {
        "status": "PASS_UNPUBLISHED_SELECTED23_EVALUATION_BOUND",
        "scope": "For U in the selected E_v only, using probability-normalized S3 measure and the positive metric Frobenius norm. Does not certify operator/source/initial closure, reduced evolution, or a numerical response. Full accepted 287-space remains at factor17 until selected-space closure is verified.",
        "pointwise_bound": "|U(n)|_F <= sqrt(321/19) ||U||_L2(S3)",
        "squared_constant": str(squared_bound),
        "sharpness": "scalar feature reproducing vector at n=e1 attains squared ratio321/19",
        "field_weights": {"scalar_yy_y0_00_trace_w": list(map(str, scalar_weights)),
                          "gradient_ya_0a": list(map(str, gradient_weights))},
        "scalar_norms": list(map(str, scalar_norms)),
        "gradient_norms": list(map(str, gradient_norms)),
        "TF_tensor_norms_HZ_HF_Tplus_Tminus": list(map(str, tensor_norms)),
        "exact_Gram_entries_checked": gram_checks,
        "TT_pointwise_squared_norms": ["1/6", "1/6"],
        "block_trace_integrals": list(map(str, expected_integrals)),
        "full23_kernel_trace_integral": "23",
        "full23_kernel_trace_at_e1": str(full_trace_at_e1),
        "full23_kernel_trace_at_e3": str(full_trace_at_e3),
        "full23_kernel_trace_polynomial": polynomial_json(full_trace),
        "certificates": certificates,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inject", choices=sorted(CONTROLS))
    args = parser.parse_args()
    need(not args.output.exists(), "OUTPUT_EXISTS")
    report = run(args.inject)
    need(args.inject is None, "CONTROL_NOT_REJECTED")
    report["checker_sha256"] = sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(report["status"])


if __name__ == "__main__":
    main()
