# Uniform nonlinear completed Cauchy data

13 November 2026. This checkpoint carries the [admitted initial constraint family](../HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/README.md) through the actual canonical-to-physical, relative-harmonic and boundary-state conversions. Under the retained smooth-background and selected weighted-family premises, the complete initial Cauchy pair is C2 on the same parameter interval and satisfies

```
|epsilon| <= 10^-6402,
sup ||(V_e(0),partial_t V_e(0))||_(H14 + H13) <= 10^-2200,
sup ||partial_e (V_e(0),partial_t V_e(0))|| <= 10^4200,
sup ||partial_e^2 (V_e(0),partial_t V_e(0))|| <= B_C=10^10400.
```

The norm is the sum of the H14 field-value and H13 first-time-velocity norms in the fixed product bundle on `Sigma=[0,6ell] x S3`, `ell=10^-6`, with full sphere measure and Frobenius tensor norms. The bounds cover one bulk copy and the retained union of two identical reflected copies. The canonical potentials remain in `X=(dot H16(q,w) scalar ⊕ dot H14(q,w) one-form)^Gamma`, `Gamma=(Z2)^4`. No angular truncation or change of the original parameter is used.

These are conservative conditional mathematical bounds in the registered normalization. Physical confirmation and observational calibration remain open. A uniform nonlinear evolution tube, observer continuation and observable remainder still need separate estimates.

## The missing conversion step now quantified

The prior checkpoint admitted exact initial constraints. Earlier Cauchy-transport estimates bounded derivatives at zero. The new calculation bounds the whole selected interval, including the prescribed lapse `alpha_e=rho(1+e f)`, `f=u_profile Y`, zero initial shift, the full mixed metric time velocities and all determinant factors in converting canonical momenta.

The actual half-action convention gives `K=2 J^-1/2[gamma tau gamma-(tr_gamma tau)gamma/3]` and normal scalar momentum `Pi=J^-1/2 omega`. Here `tau=pi/rho0`, `omega=varpi/rho0`, `rho0` is fixed, and `J=det(gamma)/det(gamma0)`. Affine canonical scalar density does not make the physical normal velocity affine.

The retained boundary straightening `F_-z` is a pointwise flow in the sixteen field variables. Its fixed-state derivative certificate applies on the proved initial patch. For `T=F_-z` composed with the actual radial square root, the initial velocity has the form `V_t=DT(P)Q`. Its second parameter derivative includes

```
DT(P) Q'' + D2T(P)[P'',Q]
          + 2 D2T(P)[P',Q'] + D3T(P)[P',P',Q].
```

The third state differential is essential even though the parameter family needs only C2 regularity. It fits the available mixed derivative budget: `14+2=13+3=16<=18`. The open initial chart patch is proved directly from the whole-family input norm, without assuming the Hessian bound being established.

The raw interfaces and displayed Hessian join in this table are for one copy. The complete proof charges factors two and four before obtaining the rounded final ceilings for the retained identical pair.

| Interface | Evaluated bound | Proof |
|---|---:|---|
| Raw completed Cauchy values | `R0<=10^-5935` | [Graded composition](GRADED_CAUCHY_COMPOSITION.md) |
| Raw first/second parameter derivatives | `R1<=10^467`, `R2<=10^6700` | [Graded composition](GRADED_CAUCHY_COMPOSITION.md) |
| Sufficient raw chart ball | `R0<=10^-3228` | [Uniform boundary transport](UNIFORM_BOUNDARY_CAUCHY_TRANSPORT.md) |
| Uniform completed Hessian join | `10^3661 R2+10^6884 R1^2` | [Uniform boundary transport](UNIFORM_BOUNDARY_CAUCHY_TRANSPORT.md) |
| Same-interval result and reflected copies | `B_C=10^10400` | [Complete proof](UNIFORM_CAUCHY_DATA.md) |

## Review, replay and next step

[AUDIT.md](AUDIT.md) records independent computational-agent mathematical review and exact proof promotion. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) describes immutable sources, closed-capsule authentication, exact algebra and arithmetic, normal/optimized execution and refreshed fault controls. These finite checks support the written analytic proof; they are not proof-assistant verification or external peer review.

[LITERATURE.md](LITERATURE.md) distinguishes pointwise field composition from varying-coordinate pullback and qualifies the scope of a relevant 2026 vacuum boundary result. [CONTINUATION.md](CONTINUATION.md) identifies the next uniform quasilinear coefficient and energy estimates. This working checkpoint creates no new Zenodo edition. Actual merge, CI and live-page receipts belong in the external handoff after publication.
