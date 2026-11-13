# Complete proper second scalar, intrinsic curvature and area at zero

This continuation follows merged PR #137. It connects the evaluated first
and second wave jets to the actual anchored proper-time observation of the
registered Einstein–scalar branch. It supplies the previously missing full
second intrinsic-curvature input and transfers it to proper area and the
normalized observable denominator at amplitude zero.

Keep `ell=10^-6`, original proper clock `q=-tau/rho_b`, original physical
anchor and transported round labels, and `0<=q<=Q=25ell/16`. Define

```
B1=10^4000000 exp(10^152017),
B2=10^10200000 exp(10^152019),
Phi2=partial_e^2 Phi_e^obs|0,
R2=partial_e^2 [Ric(h_e)(partial_q X_e,partial_q X_e)]|0,
j2=partial_e^2 j_e|0,
D2=partial_e^2 D_e|0, A2=partial_e^2 A_e|0.
```

Here `h_e` is the full physical shell metric, `j_e` the area density
relative to the background proper cut, `D_e` the retained fixed-label
Y-squared weighted denominator, and `A_e` the total proper-cut area.
Subscript two means the actual second derivative, twice the quadratic
Taylor coefficient. Under the retained analytic and geometric premises,

```
sup |Phi2| <10^7 B2+375000 B1^2 <10^8 B2,
sup |R2| <4*10^25 B2,
sup |j2|, sup |D2/D0|, sup |A2/A0| <10^14 B2.
```

These are conservative model-normalization ceilings, with no SI energy
calibration or implication of smallness or a sign of area change.

## The completed joins

[REVERSE_TIME.md](REVERSE_TIME.md) proves the actual same-state negative-time
join for the first mixed-X14 and second mixed-X13 wave bounds. It includes
the full tensor reflection and nonzero harmonic initial velocity; it does
not identify the original state's past values with its future values.
The old exact-family binding's `3ell/2` slab is insufficient here. The
separate conditional prescribed-past-slab application on `T=13ell/8`,
already stated in NONLINEAR_BRIDGE, and all its analytic premises remain
explicit conditions. The background path has margin `T-Q=ell/16`.

[PROPER_SCALAR.md](PROPER_SCALAR.md) evaluates physical shell traces and
the first observer variation in the same harmonic gauge as the second
field. The constant shell background scalar gives

```
Phi2=varphi2(X0)+2 dvarphi1(X0)[xi1].
```

[SECOND_CURVATURE_AREA.md](SECOND_CURVATURE_AREA.md) subtracts the background
Ricci relation, using `E_e=Ric(h_e)-3h_e/rho_b^2`, so `E0=0`.
With the covariant first observer velocity `eta1=nabla_e V_e|0`, it proves

```
R2=E2(V0,V0)+2(nabla_xi1 E1)(V0,V0)+4E1(eta1,V0).
```

Second observer position and velocity drop out of this identity. The
complete inverse/connection/Ricci jets, positive tensor-norm conversion,
and first observer terms are charged. The retained Raychaudhuri kernel
then supplies the full area and denominator bounds. Its curvature lower
limit remains the full family's `ell/4`, rather than the first tangent's
larger `ell` quiet interval.

## Verification and scope

[REPRODUCE.md](REPRODUCE.md) describes 119 named checks: 51 source-size/hash
checks for seventeen immutable snapshots, nine exact algebra checks
covering 119 components, and 59 rational ledger checks. The algebra
compares formulas with direct polynomial/rational differentiation,
including a curved four-dimensional Lorentzian metric with offblock
perturbations and a separate conformal-scaling control. Thirteen rejection
controls include five authentication faults and eight mathematical
coefficient faults under refreshed manifests. Normal and optimized Python
reports must agree byte for byte.

The internal independent proof review is recorded in
[INTERNAL_REVIEW.md](INTERNAL_REVIEW.md). Finite checks support the written
arguments; they do not establish their analytic existence premises or
constitute external peer review. [LITERATURE.md](LITERATURE.md) compares six
primary sources, including inspected 2025–2026 methods and their limits.

The actual family remains C2 into X10. Stronger individual jets at zero do
not make the whole parameter curve C2 into X13. Its longer past-slab
application and a qualitative common observer neighborhood remain
conditional and unquantified. No numerical nonzero amplitude interval,
uniform second derivative, finite-amplitude remainder, third amplitude
derivative, observational confirmation, or cosmological-origin result
follows. The selected odd scalar moment already has second derivative
zero under its retained symmetry premises; that is not the new output.
[CONTINUATION.md](CONTINUATION.md) states the remaining campaign.
