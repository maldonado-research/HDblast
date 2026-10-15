# Quantitative second initial data and local harmonic initialization

15 October 2026. This conditional checkpoint bounds the explicit second
initial Cauchy coefficient of the registered aligned quadrupole packet and
completes its local relative-harmonic initial gauge through parameter order
two. It does not evolve the second perturbation or establish an exact
nonlinear family.

The physical initial data were constructed in the preceding
[explicit initial constraint checkpoint](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_EXPLICIT_SECOND_INITIAL_JET_20261015).
They cancel the Hamiltonian and all four momentum sources. They are smooth,
flat at the packet endpoints, and supported in the original annulus
`ell <= z <= 2ell`, `ell=1e-6`, on both identical copies. The shell collar
initially remains the background. The potential, tension and retained
background branch are unchanged.

## Quantitative result

These strict upper bounds hold continuously over the annulus and the whole
unit three-sphere. The fields are second **parameter derivatives**, entering
the expansion as `epsilon^2 X2/2`. Bounds concern absolute individual
components in the adapted background orthonormal spatial frame, including
the derivative index. They are coarse bounds in the registered model units,
not actual response sizes, lower bounds, or physical forces.

| Physical second initial field | Component bound | Full first spatial covariant derivative bound |
|---|---:|---:|
| Spatial metric `gamma2` | `1.06e24` | `4.25e29` |
| Covariant extrinsic curvature `K2` | `7.96e7` | `1.96e13` |
| Scalar `phi2` | `5.13e33` | `3.06e39` |
| Unit-normal scalar velocity `Pi2` | `2.00e11` | `5.02e16` |

`K=+Lie_n(gamma)/2` and `Pi=n(phi)` are physical Cauchy variables, not
momentum densities. Full spatial derivatives include the angular
projector's radial-angular connection term. No Sobolev energy, nonlinear
amplitude interval, or small physical response is inferred from this table.

Two exact reductions improve the interval calculation before evaluation:
the first-order constraint removes explicit `u_y,u_yy,a_y` from the
Hamiltonian source, and the actual velocity relation cancels the large
`v v_y` terms in the scalar momentum correction. They do not change the
selected second data. An independent implementation uses unnormalized
derivatives and a different outward rounding grid to verify all 17 rounded
component caps in the detailed report.

## Initial gauge completion

Initial lapse and shift **values** are chosen as
`N=rho(1+epsilon U)`, `beta=0`, where `U=uY`. The harmonic condition fixes
their remaining time derivatives without changing the physical Cauchy
data. The explicit first and second coefficients include every quadratic
metric/lapse conversion term. Zero second lapse value does not give zero
second `g_tt` value; the mixed velocity also contains
`(dot g_ty)_2=dot beta2_y - 8 U dot beta1_y`.

The resulting full initial spacetime metric/scalar first-jet bounds, in
the fixed adapted background spacetime frame and its background covariant
derivative, are

| Second spacetime coefficient | Component bound | Full first spacetime covariant derivative bound at `t=0` |
|---|---:|---:|
| Metric `g2` | `1.06e24` | `6.40e29` |
| Scalar `phi2` | `5.13e33` | `3.06e39` |

This is a statement on the initial slice. It does not bound the spatial C1
norm of all gauge **time velocities**, or any subsequent solution. Vanishing
initial gauge coefficients and the physical constraints imply vanishing
normal gauge coefficients through order two only **conditional on the
reduced acceleration equations**. Those accelerations, the first gauge
wave, transformed bulk/shell sources and the coupled evolution remain open.
The earlier longitudinal forcing cannot be reused unchanged after harmonic
conversion.

## Derivations and independent checks

- [Physical initial bounds and source cancellations](producer/QUANTITATIVE_SECOND_INITIAL_JET.md).
- [Actual static coefficient jets](background/ACTUAL_STATIC_JET_DERIVATION.md).
- [Global bump derivatives through order six](bump/BUMP_DERIVATIVES_AND_CAUCHY_NORMS.md); order four is sufficient here.
- [Explicit local initial harmonic completion](gauge/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md).
- [Initial spacetime one-jet bounds](frame/INITIAL_SPACETIME_ONE_JET.md).
- [Independent norm review](independent-norm-review/INDEPENDENT_NORM_REVIEW.md).
- [Primary literature comparison and scope](literature/HARMONIC_INITIALIZATION_PRIMARY_COMPARISON.md).

The direct five-dimensional Christoffel check is distinct from the all-point
formal ADM expansion in `review/`. A separate all-point Laurent calculation
checks the source cancellations and the sphere bounds. The interval
certificate covers every point; generic point-jet tests supplement the
analytic proofs and do not substitute for them. Reviews are independent
internal AI reviews, not external specialist peer review or a priority claim.
The general initialization method is established in the cited literature;
this is an explicit calculation and certificate for the registered model.

## Portable reproduction

Only Python's standard library is needed. `MANIFEST.json` inventories every
source file. `REPLAY_ROSTER.json` registers the exact normal/optimized jobs,
intentional-error controls and complete expected reports. Extract the
published ZIP, use the external manifest hash pinned in
`.github/workflows/hdblast_second_initial_norms_gauge.yml`, and run

```sh
python -I -B reproduce.py --manifest-sha256 EXTERNAL_MANIFEST_SHA256 --output-dir /tmp/hdblast-second-initial-replay
```

The output directory must be new and outside the immutable source. The
runner authenticates every input, compares complete reports and verifies
the source again after all processes. The dedicated workflow separately
authenticates the ZIP and replay runner. Local source/ZIP freeze and fresh
replay evidence are published under
`hdblast/continuation/2026.10.15-second-initial-norms-gauge`.

The full second-order initial-boundary solution, exact nonlinear family,
remainder bounds and an observational test remain unfinished. No conclusion
about a higher-dimensional origin of the Big Bang follows. No new Zenodo
version or journal submission is made by this checkpoint.
