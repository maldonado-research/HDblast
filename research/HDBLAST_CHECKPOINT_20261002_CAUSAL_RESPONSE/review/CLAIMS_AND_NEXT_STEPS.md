# Scientific interpretation of the registered causal-response calibration

This is a post-run, read-only interpretation of the existing outputs. It adds
no source sample, mode evolution, quadrature grid, cutoff choice, or acceptance
threshold. `REVIEW_INPUTS.json` pins the inspected results and implementations.
The recorded public freeze is `a01d17e015f070be2ad2018745e877148fe2f4e4`;
the first primary evaluation is recorded at 2026-10-02 15:04:14.556581 UTC.

The bounded result is a successful calibration of one homogeneous retarded
scalar-variance response on a prescribed de Sitter geometry. It uses the
minimally coupled scalar at x0=r=2H^2, the fixed incoming BD state, and the
inherited positive-reference subtraction. The independent forced-mode
implementation and the finite-cutoff memory implementation agree, while the
removed-cutoff time formula agrees with the finite-cutoff results within the
separately computed analytic momentum-tail budgets.

## Evidence and practical precision

The normal and optimized validators both pass the same 36 source/time/cutoff
comparisons and reject 17 wrong-formula controls. The frozen exact verifier
passes 46 identities and detects 14 algebraic mutations in both Python modes.
These counts describe different checks and should be reported separately.

All figures below concern the registered normalized response
`y=a(eta)^2 delta Q/epsilon`, epsilon=10^-4:

| Existing diagnostic | Largest recorded value |
| --- | ---: |
| Independent modes minus primary result at the same K | 6.939e-17 |
| Fine minus coarse mode result | 1.284e-16 |
| K=256 result minus removed-cutoff memory result | 1.145e-7 |
| Analytic omitted-tail upper bound at K=256 | 2.341e-6 |
| Linear Wronskian residual divided by epsilon | 9.758e-14 |

The first two differences are much smaller than the frozen gates. They are
observed numerical agreement, not a proof of sixteen-digit accuracy. The
QUADPACK errors and mode refinement differences remain numerical estimates;
only the omitted ultraviolet tail has the stated constructive enclosure.
The tail bound is intentionally conservative, especially after the pulse.
The K=256 tail gate was 5e-6 in y. A physical variance bound is obtained by
multiplying a bound in y by epsilon/a^2.

The following continuum values are rounded from the existing primary file:

| eta | positive B: y | signed uB: y |
| --- | ---: | ---: |
| -5.5 | 0 | 0 |
| -4.5 | 0.01275119 | -0.004495970 |
| -4 | 0.0006242223 | 0.007642786 |
| -3.5 | -0.01523235 | 0.005789770 |
| -2.5 | -0.01103989 | -0.001274269 |
| -1.5 | -0.006277829 | -0.0004090006 |

Exact zero before support tests causality. The nonzero post-pulse values
test memory when the local contact vanishes. Negative post-pulse response
for B follows directly from the retarded kernel. For this particular uB
pulse its sign is also analytically negative, by pairing u and -u; this
is not a positivity theorem for arbitrary signed sources.

The state diagnostics distinguish a compact occupation perturbation and an
initial Bogoliubov perturbation from the registered vacuum source response.
They help detect an accidentally changed initial state. They do not select
a cosmological state or validate every possible excited or squeezed state.
The independent solver evolves the forced linearized equation; epsilon is
its forcing normalization. No nonlinear finite-amplitude convergence or
central-difference test was performed.

## What the result does and does not identify

The negative stationary susceptibility
`Q_x=-(2 gamma_E+ln 2)/(16 pi^2)=-0.011699927596382968` checks the normalization
and finite contact in a separate past-infinite stationary source limit. It
does not approximate the finite compact-pulse history. This distinction is
visible in the existing data: the positive pulse has positive y at its center,
whereas the instantaneous stationary replacement would give a negative value.

The scalar response can be translated into the quadratic mass-law current
only with its additional chain-rule contact:

    delta j=(br/2)delta Q+(b^2 r Q0/4)delta phi.

The renormalization contact in delta Q and this mass-law contact are distinct.
The complete translated response carries the loop/source amplitude gamma
once. This algebraic translation does not supply a stress response.

The calibrated plane-wave kernel is exact at x0/H^2=2. The previously screened
gamma=0.01 stationary root has a shifted mass and curvature and is not exactly
this reference. A causal analysis of that root requires its actual propagator
or a separately justified expansion with a controlled mismatch error. Small
gamma alone does not turn the present reference kernel into its exact Hessian.

No shell or bulk equation was evolved here. There is no completed coupled
initial-value problem, stability eigenvalue, relaxation result, physical
particle-energy yield, radiation transfer, heating or thermalization result.
The canonical spectral-work identities remain analytic companions; the run
did not integrate a spectral energy or physical stress ledger. The methods
and identities belong to established quantum field theory; no novelty or
discovery claim follows from this calibration.

## The next bounded task

Derive and separately register the matched linear response of rho, p and the
scalar current on this same fixed geometry and fixed state. The inherited
general-FRW common-action prescription already fixes the positive-reference
fourth-order stress subtraction and second-order variance subtraction. It
therefore provides a concrete next calculation without inventing new finite
counterterms from a Ward fit.

Compute rho and p independently from their minimally coupled mode operators,
including the explicit mass variation and the complete local subtraction
variations. Compare them with the independent trace-plus-Ward reconstruction
described in `stress-audit/PROSPECTIVE_MATCHED_STRESS_RESPONSE.md`. A stress
defined by solving the Ward identity cannot itself provide an independent
Ward validation. At finite comoving K, use the same finite-K reference
variance and trace remainder; inserting continuum contacts into that test
creates a spurious residual.

Before execution, freeze derivative evaluation, stress-specific ultraviolet
tail bounds, integration settings, finite-cutoff comparisons, independent
Ward/trace residuals, wrong-contact controls and a finite stopping budget.
The present Q tail bound does not automatically control differentiated Q,
pressure, or energy. Complete this fixed-geometry stress/current calibration
before extending to metric response or to a coupled shell problem.
