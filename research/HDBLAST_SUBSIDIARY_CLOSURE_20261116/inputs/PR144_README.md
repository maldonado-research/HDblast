# Reduced nonlinear evolution, full time continuity and causal support

14 November 2026 UTC. Conditional mathematical checkpoint continuing the
[uniform smooth Picard bounds](../HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/README.md)
at immutable base `9b74d57b6caef415775abc069838a4bcea303831`.

The controlled iteration now converges to an actual unique solution of
the normalized reduced sixteen-field mixed system. A separate weak-energy
and boundary argument proves strong time continuity of the entire X14
jet table. Moving-domain energy then bounds the exact solution's radial
support, keeping it separated from the artificial outer wall.

These are reduced-system results under the retained background,
coefficient-chart, selected initial-data and smooth frozen-linear
construction premises. Identifying the reduced solution with the full
unreduced Einstein–scalar/Israel problem requires subsidiary and boundary
propagation work. Physical confirmation remains open.

## Result in the registered norms

Use the full positive product bundle, all ordered covariant tensor slots,
full Frobenius metric and full sphere measure on
`Sigma=[0,6ell] x S3`, `ell=10^-6`, `T=13ell/8`. The mixed norm is

```
X_m(U;s)=sum_(j=0)^m ||partial_s^j U(s)||H^(m-j).
```

For the prescribed `Gamma=(Z2)^4` invariant finite H14/H13 Cauchy pair
`(p,q)`, supported in the stated compact-support closure on
`A0=[ell/2,5ell/2] x S3`, let `d=||p||H14+||q||H13` and require
`d<=Dbar=10^4200 a_iter`, with `a_iter=10^(-10^500209)`. The original
canonical parameter curve meets this condition when
`abs(epsilon)<=a_iter`. The constructed solution has

```
partial_s^j U in C([0,T];H^(14-j)),       j=0,...,14,
X14(U;s)<=C0 d exp(2Hs)<2*10^-15<R/4,
C0=10^7003210, H=10^500212, R=10^-12.
```

All fifteen initial jets are attained, the meaningful mixed boundary
traces hold, and the solution is independent of supported smooth
approximants and the auxiliary starting field. The original lapse,
fixed canonical density, half-action normalization, mixed velocities,
and same-state past reflection are retained. No angular cutoff is imposed
within the prescribed Gamma-invariant class. Auxiliary rotation and energy
estimates hold on all sections. The one-copy and retained identical-two-copy statements use
global input norms without a new copy factor.

The positive nested parameter ceiling is an extremely conservative
sufficient proof threshold, without a measured cosmological amplitude
interpretation. The result does not assert strong C X14 convergence of
the original approximation sequence or top-order initial-data continuity.
It does give evaluated low-energy initial-data stability and qualitative
continuity of the original parameter curve into C X13.

## Evaluated convergence and causal bounds

The new difference source has `L_delta=4*10^500000`. With
`G_en=16*10^500000`, `beta=12*10^500000`, the full-slab tail is

```
sup_s X1(U-U^N;s)
 <=(R/2) exp((G_en+beta)T) (beta T)^N/N!.
```

No short-time condition `beta*T<1` or further parameter restriction is
used. An explicit sufficient index `N_star=10^499997` gives the further
bound `R*2^(-N/2)` for `N>=N_star`. This enormous index is a theoretical
tail certificate, not a practical numerical iteration count.

For the exact initial datum supported in `[ell/2,5ell/2]`, the radial
coordinate speed is bounded by `20/19`. At the slab endpoint the support
is contained below `80ell/19`, leaving an outer margin `34ell/19`.
The perturbation remains zero at the shell before `19ell/40`; this
is a bound on possible arrival, without an arrival guarantee or a
proper-clock interpretation.

![Conditional reduced-system support envelope](figures/causal_support.svg)

The shading bounds possible support on the proved coordinate slab. It
contains no wave-amplitude calculation or observational calibration.

## Proofs and remaining scope

- [Evaluated low-energy source and factorial convergence](LOW_ENERGY_PICARD_CONVERGENCE.md).
- [Cauchy-pair compactness and exact lower mixed recovery](LOWER_MIXED_LIMIT.md).
- [Finite-data reduced existence, uniqueness and stability](NONLINEAR_LIMIT_AND_UNIQUENESS.md).
- [Closed-manifold weak linear energy regularity](WEAK_LINEAR_ENERGY_REGULARITY.md).
- [Commuted weak boundary rows, fourfold double and top time continuity](TOP_ORDER_TIME_CONTINUITY.md).
- [Integrated theorem and moving-domain causal support](REDUCED_EVOLUTION_AND_CAUSAL_SUPPORT.md).

The extra weak-energy/PDE argument establishes individual top time
continuity without an X15 estimate. Strong top approximation convergence,
top-norm data dependence, X15 tame/smoothing comparisons, evolved C2
parameters, subsidiary/Israel propagation, observer transport and a
finite nonlinear observable remainder remain separate tasks.

[Reproduction](REPRODUCIBILITY.md), [internal review/promotion scope](AUDIT.md),
[primary method sources](LITERATURE.md), and [next calculations](CONTINUATION.md)
state the verification limits. No empirical higher-dimensional origin
or new Zenodo edition is claimed by this working GitHub checkpoint.
