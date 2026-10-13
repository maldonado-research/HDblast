# Independent tensor operator-domain and exact-algebra review

13 October 2026. Review by the scalar-domain agent, independent of the author
of `TENSOR_PARTNER_GAP.md` and `verify_tensor_partner.py`. This is internal
mathematical review, not external peer review or a formal proof certificate.

**Verdict:** no blocking mathematical defect identified in the stated
conditional tensor theorem. The first-order domains, partner boundary signs,
full-spectrum comparison, kernels, and essential threshold are consistent.
The exact checker passes with and without Python optimization.

## Versions reviewed

| File | SHA-256 |
|---|---|
| `TENSOR_PARTNER_GAP.md` | `efd9d9293636864b419e46f7c333149c7847ffdc14e969d75e6a3c1e87256cad` |
| `verify_tensor_partner.py` | `bbf6cec37b1a38c4af441c59617fe647b9e44867b6350ddd5b3baf03ab673d84` |
| `TENSOR_PARTNER_CHECKS.json` | `f0515c6ffd02ef09594f7deb23bf3227b36349c5a6ba68537ebcf1619daa9f83` |
| `TENSOR_PARTNER_CHECKS_OPTIMIZED.json` | `f0515c6ffd02ef09594f7deb23bf3227b36349c5a6ba68537ebcf1619daa9f83` |

## Independent argument check

On `(-infinity,0]`, let `D=partial_z-S`, with real bounded `S` and bounded
`S'`. Integration by parts gives the adjoint expression `D*=-partial_z-S`.
Taking `Dom D=H^1` makes `Dom D*` the zero-trace `H^1` space; reversing these
domains reverses that endpoint assignment. No cone boundary term survives
for the half-line Sobolev domains. These are closed densely defined operators.

For the even assignment, `D*D` has `D psi(0)=0`, and `DD*` has Dirichlet
trace. For the odd assignment, `D*D` has Dirichlet trace, and `DD*` has
`chi'(0)+S_b chi(0)=0`. Expanding the latter partner form gives

    ||D* chi||^2 = integral (|chi'|^2+(S^2-S')|chi|^2)
                  + S_b |chi(0)|^2.

The positive sign is correct for the specified coordinate orientation. Thus
`U<=0` and `S_b>0` imply `DD*>=9/4` in both assignments, using

    S^2-S' = 9/4 + rho^2 (9 p^2/16 - U/8).

Polar decomposition compares the entire nonzero spectral subspaces of
`D*D` and `DD*`; this is not restricted to eigenfunctions. Consequently no
spectrum lies in `(0,9/4)`, and factorization excludes negative spectrum.
The only even kernel is `rho^(3/2)`, whose squared norm is `integral rho^2 dy`
and is finite. Its nonzero shell trace excludes it from the odd domain.
The formal partner kernel `rho^(-3/2)` has divergent cone norm.

The absence of an odd kernel also excludes zero from the odd spectrum:
a possible spectral point at zero would be isolated by the proved open gap,
and an isolated spectral point of a self-adjoint operator is an eigenvalue.
For the even operator, zero is instead a simple isolated eigenvalue.

The smooth regular-cone bounds make the potential differences from `9/4`
bounded and decaying on the half-line. Local Rellich compactness plus a
small-norm tail establishes relative compactness of the potential
perturbation. The regular endpoint change has finite-rank resolvent effect.
Equivalently, polar decomposition transfers the partner essential spectrum;
the even zero eigenspace is finite dimensional. Thus the essential spectrum
is `[9/4,infinity)` for both original parities.

## Exact-check replay and separate arithmetic

Independently reran the author's checker with normal Python and `python -O`,
writing receipts outside the source directory. Each run passed eight
identities/sign checks and rejected four deliberately incorrect formulas.
The receipts are byte-identical and match the pinned receipt hash above.
The checker uses explicit failures, not optimization-removable assertions.

A separate exact `fractions.Fraction` calculation reconstructed the displayed
coefficient bound on `|eta|<=1/100` and obtained

    U_upper = -997874976899/13500000000000 < 0.

It also independently checked the coefficient reductions
`(3/4)(1/12)+1/2=9/16` and `3/4+3/2=9/4`. The polynomial itself agrees by
expansion with `U=(W')^2/2-2 W^2/3`, `W=1/3+eta^2+eta^3/3`.

These computations check exact algebra and a polynomial interval inequality.
They do not mechanically verify the functional-analytic argument or physical
mode reduction.

## Scope

The mathematical conclusion is conditional on the stated pure transverse-
traceless radial equation, parity junctions, normalizability choice, and
certified background properties. This review did not independently rederive
the full gravitational second variation or doubled-shell tensor junction.
The physical reduction and exceptional-harmonic overlap remain the subject
of the separate physical audit. No scalar kinetic-positivity, vector-sector,
nonlinear-stability, cosmological, or uniform positive physical-mass theorem
is inferred here. The distinction between the dimensionless threshold
`9/4` and the physical scale `H_b` is maintained in the reviewed note.
