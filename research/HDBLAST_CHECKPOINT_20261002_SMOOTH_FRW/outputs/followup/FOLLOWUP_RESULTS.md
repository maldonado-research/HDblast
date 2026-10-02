# Separate K384 follow-up: executed result

The separately registered curved follow-up **passes all 13 acceptance groups**.
Its K192-to-K384 pressure change is 0.005086914388911179, below the unchanged
allowance 0.011121896652637492. This is a change of 0.229101% of the larger
maximum pressure signal. The original K192 matrix remains **FAIL**; its
pressure change exceeded its original cutoff allowance.

The new run preserves the same A=0.2,r=4,T=1 prescribed smooth background,
exact incoming state, complete stress/current subtraction, finite-action
convention and physical solver. Only the separately registered momentum
range was extended. Four new runs tested eight versus twelve quadrature
nodes, phase_step=0.1 versus 0.05, and an exact-static K384 negative control.
Both new K192 subsets match their corresponding preserved original curves.

## Numerical acceptance evidence

| Maximum-over-time comparison | rho | p | Q |
|---|---:|---:|---:|
| K192-to-K384 change | 1.67988e-5 | 5.08691e-3 | 1.52253e-5 |
| Same-cutoff phase-step refinement | 1.42529e-6 | 5.20675e-7 | 1.63368e-11 |
| Same-cutoff quadrature refinement | 6.63911e-7 | 2.38023e-7 | 1.02147e-11 |
| Exact-static K384 residual | 1.74487e-6 | 5.99853e-7 | 1.77074e-11 |

The same ordinary observable allowance is 2e-5 plus 0.5% of the larger global
signal maximum. It is not a uniform pointwise relative-error guarantee near
zero crossings. For the cutoff comparison, the normalized rho,p,Q changes
are respectively 0.0119716%,0.229101%,0.0239631%. No asymptotic power was fitted
or used to accept the result, and no K>384 mode was executed.

The selected twelve-node K384 run has raw normalized Wronskian error
8.77076e-14. Its maximum integrated paired-work ledger residual is 5.29197e-6,
well below the registered allowance 1.56279e-3; the independently accumulated
ODE work and the 401/801 time-ledger comparison also pass. Its local sampled
exchange residual is 6.18692e-4 versus allowance 1.01809e-2.

The independent finite-K trace moment check gives 1.26406e-11 with direct
physical-equation Q derivatives and 1.45881e-5 with sampled Q derivatives.
The sampled 401-node trace residual is 2.29697e-4. These pressure checks pass
alongside the independent cutoff, quadrature and solver comparisons.

Cancellation of large vacuum terms is measurable: the K384 static rho floor
is approximately 1.7e-6, and the future full-rho/spectral-energy difference
is 2.37773e-6. Increasing the cutoff reduced the pressure tail while increasing
the floating-point floor relative to K192. These observed effects are retained
in the diagnostics rather than removed by enforcing a Wronskian normalization.

## Coherent state and paired work

The selected future full-minus-occupation pressure has maximum magnitude
0.0505475; the corresponding Q difference is 0.0234771. Full coherent pressure
and current were retained. Omitting source work raises the integrated ledger
residual to 0.211969, compared with 5.29197e-6 when the paired x'Q/2 term is
included. These results support the implemented exchange balance on the
prescribed background, not an energy-consumption calculation in a coupled
background evolution.

All sampled arrays are finite. At the x=0 crossing eta=0.5, the selected K384
values are rho=-0.02371166027,p=0.16335385841,Q=0.02655551722. At eta=1.25,
rho=0.13695081219,p=0.03582304839,Q=0.00510435148. The exact future occupation
energy is approximately 0.13695272808; its difference from direct subtracted
rho is covered by the recorded cancellation diagnostics.

## Registration and scope

The follow-up was committed before execution as
dd1f00fc62a5b2d257f47b8490b79de1514bd2f3. Its protocol SHA256 is
b067afda194895255081553e7bd9a3833b781f5277d3011b619829b4b3033d51.
The wrapper SHA256 is
334e0f347df0782765edbbaba891f79f8d8bc48ce9d58275804bd843e8e4c430,
and it verifies the unchanged repaired core SHA256
32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a.
The compact baseline carries exact original summary/curve hashes and selected
column equality against the retained full original arrays. Every threshold
and formula remained unchanged after registration.

This directory contains final diagnostics and every integrated curve from the
four new runs. Its curation manifest hashes these files and eight retained full
raw artifacts. Numeric NPZ files use no object arrays and are loaded with
allow_pickle=False. Plot curves come directly from the executed arrays.

A pass establishes finite-cutoff numerical agreement under this prescribed
smooth-background experiment and the declared scheme, at the registered
tolerances. It does not prove the infinite-cutoff limit, change the earlier
FAIL result, alter the archived shell, or demonstrate backreaction, decay,
thermal radiation, sustained expansion, or a new cosmological solution.
