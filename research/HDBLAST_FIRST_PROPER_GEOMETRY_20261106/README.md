# Evaluated first proper-cut geometry

6 November 2026 UTC (5 November in Los Angeles). This checkpoint evaluates
actual first variations of the selected shell's proper-time cut metric and
area density. The amplitude is the fixed registered affine parameter e;
derivatives are evaluated at e=0, with original proper time and transported
initial labels held fixed.

On `0<=q<=25ell/16`, `q=-tau/rho_b`, `ell=10^-6`, the bounds are

```
sup ||G'_0||op < 2.4*10^-8,
sup |j'_0|    < 3.5*10^-8,
G_e=k_0^(-1/2) k_e k_0^(-1/2),   j_e=sqrt(det G_e).
```

These are derivatives per unit registered amplitude, not finite-amplitude
percentage changes. They use the retained rigorous actual master-trace
certificate and same-state past/future transfer, followed by a new
first-variation calculation. The cut metric includes the proper-clock
shift and angular observer transport. Jointly transporting fields,
physical anchor and observers removes the tangential harmonic gauge
terms from these physical pulled-back coefficients. Longitudinal
coordinate clock/angle bounds are recorded separately and cannot be used
as harmonic-coordinate U1 or uniform nonlinear P1/P2.

The derivation also checks that the transported-template denominator has
zero first variation, consistent with the earlier quarter-turn symmetry.
Its first-tangent no-response interval `q<=ell` must not replace the full
nonlinear family's shorter protected interval `q<=ell/4`.

This supplies evaluated input `P0=sup||G'_0||<2.4*10^-8` to the nonlinear
campaign. If a uniform second cut-metric bound `Q(a)=sup||G''_e||` is later
validated on the same admitted parameter/domain box, then

```
sup ||G'_e|| <= P0+a Q(a),
sup ||G_e-Id|| <= a P0+a^2 Q(a)/2.
```

Q(a), the actual admissible radius, full nonlinear scalar/metric/embedding
and observer estimates remain unevaluated. Small first variations alone
do not establish finite-amplitude spacelikeness, noncrossing, a remainder
bound, or the interacting central zero. A normalized abstract metric
counterexample in the proof explains this logical gap; it is not a
solution of the selected Einstein–scalar system.

Read [the continuous derivation](FIRST_CUT_BOUND.md),
[immutable source pins](SOURCE_PINS.json), [method references](REFERENCES.md),
and [reproduction instructions](REPRODUCE.md). The arithmetic checker
consumes the pinned historical actual-field reports; it does not rerun
their full producer or certify the retained PDE/reconstruction premises.
Independent internal review is distinct from external peer review.

This is a focused GitHub working checkpoint, without a new Zenodo version
or journal submission. No higher-dimensional origin of the Big Bang or
observational confirmation of HDBLAST follows.
