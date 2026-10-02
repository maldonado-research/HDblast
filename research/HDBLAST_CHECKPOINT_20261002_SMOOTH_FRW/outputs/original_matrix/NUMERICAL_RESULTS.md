# Registered smooth FRW control: executed numerical result

The registered matrix **fails overall**. The independent flat control (A=0)
passes every final acceptance gate. The curved control (A=0.2) fails its final
pressure cutoff-convergence gate; its other final gates pass. The exact-static
negative control also passes. No tolerance was changed after viewing results,
and the registered maximum cutoff K=192 was respected.

The executed physical problem is r=4,T=1 with a=exp(A B), x=4(2B-1)^2, the
specified C-infinity step, and an exact static incoming vacuum. This is a new
prescribed-background control. It does not alter the archived shell or solve
quantum backreaction, decay, or thermalization.

## Main acceptance results

All differences are maxima over eta in [0,1.25]. Both amplitudes invoked the
prospectively registered K=192 extension because pressure failed at K=96.
The final selected curves use twelve Gauss-Legendre nodes per width-two panel,
DOP853 rtol=2e-13,atol=2e-15, and phase_step=0.1. Separate eight-node,
phase_step=0.05 integrations tested solver sensitivity at K=192.

| Quantity | Flat A=0 | Curved A=0.2 |
|---|---:|---:|
| Final matrix status | PASS | FAIL: pressure cutoff |
| K96-to-K192 max rho change | 9.53325e-6 | 7.39685e-5 |
| K96-to-K192 max p change | 8.29198e-4 | 2.26243e-2 |
| Pressure change allowance | 3.03765e-3 | 1.12327e-2 |
| K96-to-K192 max Q change | 2.26815e-5 | 6.20191e-5 |
| Max raw normalized Wronskian error, selected solve | 5.11813e-14 | 3.81917e-14 |
| Integrated paired-work ledger residual | 1.82080e-7 | 1.45806e-7 |
| Sampled local exchange residual, 801 times | 6.68238e-5 | 6.73505e-6 |
| Finite-K trace residual, analytic mode derivatives | 2.31815e-13 | 8.45990e-13 |
| Finite-K trace residual, sampled Q derivatives | 1.01195e-6 | 1.48555e-5 |
| Future spectral-energy versus full-rho difference | 1.82000e-7 | 6.54801e-8 |

The curved pressure cutoff change is 1.008866% of the larger signal maximum,
exceeding the registered 0.5% plus 2e-5 allowance. Its K192 quadrature change is
1.60342e-8 and its phase-step solver change is 2.38075e-8. The cutoff discrepancy
is therefore much larger than either measured numerical integration change.
Passing exchange and trace identities does not remove this discrepancy.
The curved result is a finite-cutoff calculation with unresolved pressure
convergence at the registered tolerance.

The exact-static control gives maximum absolute residuals rho=7.92459e-9,
p=2.73647e-9,Q=1.17390e-12. At K192 the selected future spectral-energy residual
is approximately 1e-7: cancellation of large vacuum terms remains measurable
even when individual mode Wronskians are accurate to approximately 1e-14.
The finite-difference local exchange residual need not improve monotonically
with denser sampling once this cancellation floor is differentiated. Both
resolutions and the independently accumulated work are retained.

## Paired current and coherent stress

The integrated exchange ledger explicitly includes x'Q/2. Omitting this term
raises the maximum ledger residual to 0.127790 in the flat control and 0.211980
in the curved control, compared with approximately 1e-7 with the paired work.
This negative diagnostic demonstrates that the scalar current is needed in
the measured exchange balance.

In the exactly static future, full rho agrees with the exact occupation energy
within the registered allowance. Pressure and Q retain coherent contributions:
the maximum full-minus-occupation pressure is 0.147421 for A=0 and 0.0505502 for
A=0.2; the corresponding Q differences are 0.0156263 and 0.0234771. These terms
were retained throughout the calculation. Occupation-only pressure or Q would
not represent the computed state.

All mode and integrated samples are finite, including at x=0. At eta=0.5 the
selected curved finite-cutoff values are rho=-0.0237050,p=0.163381,Q=0.0265494.
They are reported as finite-cutoff samples, subject to the pressure gate above,
and carry no claim about energy consumption by a dynamical background.

## Execution record and reproducibility

The prospective protocol was committed as
983c472dfe70d2a830b8e9b538995732dc38626f, with SHA256
37fda5b30332d0079ec44db708e266788fdf5859303707bb709750fdcf59d898.
The first registered execution saved its physical modes and then failed when
NumPy-left arithmetic bypassed the custom Taylor-jet reflected operators.
The failure, original source and physical modes were preserved.

The mechanical repair added Jet.__array_priority__=1000; it did not change
formulas, background, matrix, or acceptance thresholds. It was committed as
57a97a281879f03f53a127cd6f1f00ddb3c55eb7 before the fresh run. The repaired
source SHA256 is
32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a.
A no-mode preflight verifies ndarray-left operations, shapes, finite
counterterms, the static/flat formulas, and the complete subtraction exchange
identity on selected samples. The independent Radau and symbolic reviews are
recorded in separate artifacts.

This directory contains the unchanged final summary, per-run diagnostics and
all integrated curves. Numeric NPZ files must be loaded with allow_pickle=False.
The curation manifest hashes these outputs and all 26 retained raw mode/array
artifacts, which remain locally preserved and are reproducible from the
registered scripts. Figure curves are drawn directly from executed arrays.

Any further K>192 research requires a separately registered follow-up. It must
preserve this matrix's FAIL classification and independently control solver,
quadrature and cancellation errors. The existing result does not justify
declaring the curved pressure converged or importing it into shell evolution.
