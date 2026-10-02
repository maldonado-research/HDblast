# Actual massive de Sitter quantum sources: four-point mathematical benchmark

The registered determinant calculation passed all 57 primary gates and detected
three deliberately wrong source prescriptions. The separately implemented
proper-time calculation agrees at all four points and passes all 12 registered
cross-method comparisons. These are one-loop Euclidean/Bunch–Davies vacuum
sources in the explicitly declared finite convention, with no physical HDBLAST
mass or gravitational coupling selected.

| x/r | H²/r | W/r² | rho/r² | Q/r |
|---:|---:|---:|---:|---:|
| 1 | 0.5 | -5.431202305394E-4 | 2.902429739754E-4 | 4.221715985097E-3 |
| 1 | 1 | -4.411802236985E-3 | 2.119232911866E-3 | 2.071990800425E-2 |
| 2 | 0.5 | 1.280919704858E-5 | -4.118184790862E-5 | -1.351351456310E-4 |
| 2 | 1 | 5.592079284983E-5 | 3.694001486960E-4 | 2.110857992549E-3 |

The pressure is p=−rho by de Sitter invariance. The scalar-shell source remains
J_phi=x_phi Q/2; no numerical x_phi is invented. The negative rho and Q at
(x/r,H²/r)=(2,1/2) are retained. They are finite-scheme vacuum expectation
values, not a negative occupation number or a particle-heating result.

The calculation derives both sources from the same finite action: Q=2W_x and
rho=W−(H²/2)W_(H²), with r fixed. The independent calculation instead integrates
the heat trace and its separately differentiated mass and metric kernels.
It does not import the zeta producer. The trace relation is checked only after
the metric density has been computed; it is an algebraic consistency check of
the shared action and is not independent evidence for determinant truncation.

## Numerical evidence and its limits

The table gives the largest absolute quantities across the four points, in
the same reference units as the source table. The two primary calculations
use 50 decimal digits/N=128 and 80 decimal digits/N=256, respectively.

| Observable | Primary analytic series-tail bound | Refined analytic series-tail bound | Precision/truncation change | Refined proper-time difference |
|---|---:|---:|---:|---:|
| W | 1.561E-18 | 8.360E-33 | 1.730E-19 | 1.015E-31 |
| rho | 2.332E-16 | 2.472E-30 | 2.542E-17 | 5.986E-31 |
| Q | 4.602E-16 | 4.911E-30 | 5.083E-17 | 5.441E-31 |

The explicit geometric bounds control the convergent zeta-series tail only.
They do not certify mpmath transcendental evaluation or floating-point roundoff.
The independent proper-time calculation has analytic harmonic and infrared
tail bounds, but its small-time asymptotic remainder and quadrature error are
tested by refinement without a proved total-error enclosure. Thus the agreement
is strong numerical evidence within the registered absolute/relative gates,
not an interval-certified claim of every displayed digit. No tolerance or grid
was changed after seeing the sources.

Independent centered action differences with Richardson extrapolation agree
with Q to 1.91e−19 and rho to 1.33e−21 at worst. The independently differenced
common-action pairing residual is at most 5.76e−21; the separate digamma Green
function differs from Q by at most 5.45e−31. Negative controls omit the metric
variation, reverse the scalar current, and drop the current paired to a local
F(x)R term; all are rejected at the original gates.

## Prospectivity and scope

The source baseline is 0205cc651bfb614c32e39dfe833d93d957264229. The local
registration and code hashes were frozen at 2026-10-02T06:17:06.731272Z,
before any primary source evaluation. The identical files were publicly
committed/pushed as 85aea9955b99bba911a0e66e5869c184dd86a660 before the
original primary run began at 06:18:58.957725Z and completed at 06:18:59.751705Z.
The independent reviewer ran under an earlier local registration; its run was
completed before this public commit. These distinct timing claims must be
preserved. The registered producer, validator, source inputs and run provenance
are recorded by SHA-256. No primary run failed or required a repair.

This is a bounded mathematical source benchmark. It does not establish a
physical shell parameter choice, solve either junction or a coupled radial
boundary-value problem, evolve the bulk or shell, prove quantum stability,
or calculate heating or particle production. Its four points are not an
extrapolation over the full mass/curvature domain.
