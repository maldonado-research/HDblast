# Evaluated mixed evolution estimates and smooth Picard bounds

14 November 2026 UTC. Conditional mathematical checkpoint continuing the
[initial jets and compatible seed](../HDBLAST_TEMPORAL_JETS_SEED_20261114/README.md)
at immutable base `d27cfd999aa2e57995c08b1228a73d3da7c973d1`.

The actual frozen sixteen-field equation now has evaluated high-order
energy and interior normal-recovery constants. Their join supplies a
uniform X14 bound for compatible smooth Picard iterates on the prescribed
past slab. The finite-Sobolev nonlinear limit and evolved observer results
remain open.

The norms retain the fixed positive product bundle, full tensor/Frobenius
derivative slots and full sphere measure on `Sigma=[0,6ell] x S3`,
`ell=10^-6`, `T=13ell/8`. Define

```
X_m(w;s)=sum_(j=0)^m ||partial_s^j w(s)||H^(m-j),
I14=sum_(j=0)^14 ||J_j||H^(14-j),
Y14(u;s)^2=sum_(p+length(I)<=13) E_v[partial_s^p R^I u](s).
```

The six natural rotation generators act on every tensor slot; words are
ordered. Their estimates apply to all sections, including words outside
the original `Gamma=(Z2)^4` invariant subspace. No angular cutoff or
invariant spectral gap is used. One copy and the retained identical-two-copy
Hilbert union have the same constants when the input norms already use
that global union.

For a smooth compatible predecessor with `sup_s X14(v)<=R=10^-12`
and arbitrary smooth frozen successor u,

```
Y14(u;s)<=C_Y I14+integral_0^s[K_u X14(u;r)+K_v X14(v;r)]dr,
X14(u;s)<=N Y14(u;s)+L_u X13(u;s)+L_v X13(v;s).
```

| Constant | Evaluated ceiling |
|---|---|
| `C_Y` | `10^30` |
| `K_u`, `K_v` | `10^500101`, `10^500100` |
| `N` | `10^110` |
| `L_u`, `L_v` | both `10^500208` |

The predecessor hypothesis places its coefficients inside the actual
state patch. No successor smallness is assumed in these two estimates.
The [normalized metric margins](NORMALIZED_METRIC_MARGINS.md),
[rotation and mixed-product proof](ROTATION_MIXED_PRODUCTS.md),
[energy/commutator proof](FROZEN_ENERGY_COMMUTATORS.md), and
[normal-recovery proof](INTERIOR_NORMAL_RECOVERY.md) account for every
coefficient derivative, boundary flux, word and normal stage. Normal
derivatives are recovered from the unisolated interior equation without
normally differentiating boundary conditions or assuming inverse-coefficient
jet bounds.

The [joined smooth-iteration estimate](CONDITIONAL_PICARD_ADMISSION.md)
uses `C0=10^7003210`, `H=10^500212` and the sufficient restriction
`abs(epsilon)<=a_iter=10^(-10^500209)`. This positive nested power is
strictly smaller than the preceding seed interval. It is a conservative
proof threshold in the same original parameter, without a measured
physical amplitude interpretation. Allowing a factor two for supported
smooth approximants gives `sup_(n,s) X14(U^n;s)<2*10^-15<R/4`.
The common initial jets remain J0 through J14.

These statements retain the actual smooth-background, coefficient/flow
chart, canonical initial-family and smooth frozen-linear construction
premises. Uniform smooth-iterate bounds alone do not prove convergence,
continuous full X14 completion, evolved C2 dependence or a nonlinear
observer tube. A physical preparation, calibration and observational
test remain separate work. This checkpoint supplies no empirical
higher-dimensional event or finite nonlinear observable remainder.

[Reproduction](REPRODUCIBILITY.md), [review and promotion scope](AUDIT.md),
[primary method sources](LITERATURE.md), and [next calculation](CONTINUATION.md)
record the verification boundary. No new Zenodo edition accompanies this
working GitHub checkpoint.
